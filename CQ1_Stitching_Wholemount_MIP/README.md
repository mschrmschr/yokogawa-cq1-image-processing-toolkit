# CQ1_Stitching_Wholemount_MIP

Stitches raw CQ1 (Yokogawa) wholemount tile scans into per-well/per-grid
FIJI-compatible OME-TIFF mosaics, and automatically writes a colorized
per-channel max-intensity-projection (MIP) panel PNG alongside each one, for
a quick-look overview (not another OME-TIFF — there is no separate MIP tif).

This is a cleaned-up, working port of the stitching pipeline found (in various
broken/partial states) throughout `CQ1_Processing_scripts_noFIJIcompatbility`.
The most complete prior version was `CQ1_Stitch_99%_MIP`; this folder fixes
two bugs found there (see "Fixes vs. the source scripts" below).

## Why this is "FIJI compatible"

Every OME-TIFF this pipeline writes goes through
`tifffile.imwrite(path, array, ome=True, metadata={"axes": ..., "Channel": {"Name": [...]}, ...})`
and **never** passes a hand-built `description=` string. That distinction
matters: `tifffile` only auto-generates real OME-XML from the array's actual
shape when you let it build the metadata itself via `metadata=`; passing your
own `description=` instead makes it silently fall back to a generic "shaped"
TIFF description that Bio-Formats/FIJI cannot open as a multi-channel Z-stack
hyperstack. Keep that rule in mind if you extend this code.

## Expected input layout

```
<root_dir>\
  MeasurementResult.ome.xml     CQ1's own OME-XML export
  Image\                        raw per-tile TIFFs referenced from the XML
```

`ome_metadata.py` parses the XML into a per-tile CSV (`WellID`, `GridIndex`,
`Z`, `C`, stage position, physical pixel size, and both lattice and
measured-stage pixel coordinates). It's generated automatically on first run
and cached at `meta_csv` (default: `<root_dir>\image_metadata_with_stage_and_grid.csv`).

CQ1's own OME-XML only gives generic per-channel labels (e.g. "CH1", "CH2"),
not real staining names — it doesn't know what you stained with. Set
`channel_names` in `jobs.json` (a list, one per channel, ordered by
ascending `C` plane index) — e.g. `"channel_names": ["DAPI", "GFP", "RFP"]`
— to get real names anywhere they matter. If omitted, the generic "CH1"-style
XML labels are used as-is, so for anything other than a quick test always set
this explicitly (the GUI requires it for this reason). This renames the
channels in the stitched OME-TIFF itself, and the MIP panel PNG picks up the
same names automatically (it reads them back out of the stitched file). You
can still set `mip.channel_names` separately if you want the MIP panel
labeled differently from the stitched output. On the CLI, use
`--channel_names DAPI GFP RFP`.

**All wells in one `jobs.json` job share the same `channel_names`** — there's
no per-well override. If different wells in the same dataset were stained
differently, run this tool once per staining group instead, using "Only well"
(GUI Advanced options; `only_well` in `jobs.json`/CLI) to scope each run to
just that well and setting `channel_names` accordingly each time. Output
filenames already include the WellID, so repeated runs into the same output
folder won't overwrite each other's stitched results — and when `only_well`
is set, the MIP panel PNG step is automatically scoped to that well's files
too (via an auto-added `include_keywords` filter), so a later run's
`channel_names` can't leak into an earlier well's already-generated panel PNG.

## Running

1. `pip install -r requirements.txt`
2. Edit `jobs.json` — one object per dataset. At minimum set `root_dir`.
3. `python run_jobs.py jobs.json` (or just `python run_jobs.py`, which
   auto-discovers `jobs.json` next to the script).

Each job stitches every `(WellID, GridIndex)` pair found in the metadata CSV,
then writes a colorized MIP panel PNG for whatever it just wrote. A pair for
which any listed tile file is missing is not written; the log names the files.

### Output

```
<root_dir>\stitched_output\
  W0001_A1_stitched.ome.tif        stitched wholemount, axes ZCYX
  W0001_A1_stitched_placement.json how the tiles were placed (measured method)
  ...
  mip\
    W0001_A1_stitched_MIP_panel.png   colorized per-channel MIP panel, quick-look only
    ...
```

The stitched wholemount is a standard OME-TIFF openable directly in FIJI via
Bio-Formats (`File > Import > Bio-Formats`), with correct channel names and
pixel calibration.

### Output compression (`compression` / `compression_level` in `jobs.json`,
"TIFF compression" / "Compression level" under Advanced options in the GUI)

Every stitched stack is already written with **lossless** compression by
default (`compression: "zlib"`, `compression_level: 6`) -- no pixel values
are ever changed by any of these settings, only how compactly they're
packed on disk.

Choices, benchmarked against a real 1.62 GB stitched stack (20 sampled
planes, `(3789, 3789)` uint16 each):

| `compression` | size vs. raw pixels | vs. `zlib` default | opens in FIJI/Bio-Formats |
|---|---|---|---|
| `zlib` (default) | 44.7% | -- | yes |
| `zstd` (level 19) | 40.1% | ~10% smaller | yes -- verified byte-identical read-back |
| `none` | 100% | much larger | yes |

`zstd` is the one worthwhile alternative: genuinely smaller, still fully
lossless, and confirmed (via a headless Fiji/Bio-Formats round-trip, not
just `tifffile` reading its own output) to decode back byte-for-byte
identical to the source. It's slower to write than `zlib`, though writing
is a small fraction of a stitching job's total time. `compression_level`
means something different per scheme: `zlib` uses 0-9, `zstd` uses 0-22
(higher is smaller/slower for both).

`LZW` and `LZMA` were also tried and rejected: `LZW` compresses *worse*
than `zlib` on this data (CQ1 fluorescence images are dominated by
photon shot noise, not the kind of redundancy LZW exploits), and `LZMA`
compresses best of all (~36% of raw) but **isn't readable by Bio-Formats
at all** (`EnumException: Unable to find TiffCompresssion with code:
34925`) -- confirmed by trying to open an LZMA-compressed sample in Fiji
headless before it went anywhere near a real dataset. Neither is exposed
as a choice.

The TIFF "predictor" option (horizontal differencing before compression)
is deliberately a no-op regardless of the `predictor` setting -- see the
comment on `_writer_kwargs()` in `stitcher_unified.py`. It measurably
helps on smooth synthetic gradients but measurably *hurts* on this
project's real, noise-dominated data (~2.6% larger on the same benchmark
file). Left in place as a no-op rather than removed, so a future config
change here doesn't need to relearn this.

## Stitch backends (`stitch_backend` in `jobs.json`)

- **`measured`** (default; "Tile placement: measured" in the GUI) — tile
  positions are measured from the image overlaps, see below.
- **`classic`** — straight tile paste, no blending. Fastest; visible seams if
  overlap/positioning isn't exact.
- **`seamless`** ("Tile placement: stage" in the GUI; the default of earlier
  versions) — feather-blends overlapping regions with a
  raised-cosine taper; supports optional lattice snapping, local overlap-based
  alignment refinement, and overlap gain-matching. The taper (here and in
  `tiny`) is applied only where an earlier tile has already written pixels, so
  the outer mosaic edge keeps the tiles' own values. (Before 2026-11-25 it
  faded to black there, by up to the feather width.)
- **`tiny`** — for datasets with near-zero (~0-1%) nominal tile overlap:
  snaps tile positions to the measured lattice and gain-matches each tile to
  the grid's median brightness before feather-blending a narrow seam.

`classic`, `seamless` and `tiny` place tiles from the CQ1's recorded stage/grid
positions (`placement: "stage"` or `"grid"` in `jobs.json`). On real
acquisitions those leave neighbouring tiles several pixels apart, by an amount
that differs from join to join.

### Measured placement (`tile_placement.py`)

1. For every pair of neighbouring tiles, the two facing edge strips (maximum
   projection of up to 8 planes, local mean removed) are matched by normalised
   cross-correlation, within +-60 px along the join (at most a third of the
   overlap) and +-40 px across it. Every shift in that window is scored exactly
   (Fourier transform for the sums of products, summed-area tables for the rest:
   J. P. Lewis, "Fast Normalized Cross-Correlation"). Tiles are read 8 at a
   time, which matters on network drives.
2. A join counts as measured if at least 2 channels reach a correlation of 0.5
   and all of those agree within 3 px. With a single channel it needs 0.8 and
   agreement with another join of the same direction. Overlaps that leave less
   than 32 px to compare are not measured.
3. A join that is not measured gets the typical (median) offset of the measured
   joins in the same direction. If there is none, so that a join would stay at
   the stage position, **the whole grid is stitched with the stage method**
   (`seamless`, narrow blend) instead; the log says so and the JSON records
   `"method": "stage (fallback)"`. Single-channel datasets often end up here.
4. Positions are solved from all joins together; each tile gets one position for
   all planes and channels.
5. Blending: every pixel is a raw value or a weighted mean of raw values. The
   weight passes from one tile to the other in the middle of the overlap over
   `blend_px` pixels (default 128, never more than the smallest overlap).
   Brightness differences between tiles are spread over that width, **not**
   corrected.

Next to every stitched file, `<name>_placement.json` records each join (measured
or not, offset, correlation per channel, leftover disagreement) and the final
tile positions; the log summarises it. The feather, lattice-snap, refine and
gain settings do not apply to this method (only to a grid that falls back).

### Sources (measured placement)

Used:

- J. P. Lewis, "Fast Normalized Cross-Correlation" (Industrial Light & Magic),
  <https://scribblethink.org/Work/nvisionInterface/nip.html>. The page lists
  an earlier paper: J. P. Lewis, "Fast Template Matching", *Vision Interface*,
  1995, pp. 120-123. Used for: sums of products by Fourier transform, window
  sums and sums of squares from summed-area tables. Lewis describes a fixed
  template searched in an image; here both overlap windows change with the
  shift, so both images get summed-area tables.
- F. C. Crow, "Summed-area tables for texture mapping", *Computer Graphics*
  (SIGGRAPH '84) 18(3), 207-212, <https://doi.org/10.1145/800031.808600>. The
  summed-area table.

The accept rules, the least-squares tile positions and the blend weights are
this project's own; their thresholds were set on its data and are listed above.

Related, not used:

- D. Padfield's masked normalised correlation in the Fourier domain (IEEE
  Transactions on Image Processing, 2012), as implemented for example in
  [ITK](https://docs.itk.org/projects/doxygen/en/stable/classitk_1_1MaskedFFTNormalizedCorrelationImageFilter.html).
- [MIST](https://github.com/USNISTGOV/MIST/wiki) (NIST), a stitching tool that
  uses phase correlation followed by normalised cross-correlation.

## Fixes applied vs. the source scripts

1. **Backend/kwargs mismatch.** In the source `CQ1_Stitch_99%_MIP`, the batch
   driver always builds one large `blend_args` dict and forwards all of it to
   whichever backend was selected — but `stitch_grid_classic` only accepted
   `placement`/`verbose`, and `stitch_grid_tinyoverlap` didn't accept several
   keys `run_jobs.py` always sets (e.g. `lattice_snap_px`, `refine_*`,
   `overlap_gain_match`). Only `seamless` happened to accept everything, so
   `classic`/`tiny` crashed with `TypeError` when selected. Here, `main.py`'s
   `_dispatch()` filters `blend_args` through each backend's actual signature
   (via `inspect.signature`) before calling it, so all three backends work
   through the same job driver.
2. **Dead `tile` option.** The source scripts accepted and even documented
   `tile` in `jobs.json`, but never passed it into the
   `tifffile.imwrite(...)` calls — it had no effect. Here it's wired in: a
   tiled BigTIFF is important for opening large wholemount mosaics in FIJI
   without loading the whole image into memory. `pyramids` is still accepted
   in `jobs.json` but remains a no-op (as it was in the source scripts) —
   real pyramid SubIFDs need tifffile's lower-level `TiffWriter` API plus
   explicit per-level downsampling, which the installed tifffile version's
   `imwrite()` convenience wrapper doesn't expose directly.
