# tile_placement.py - the "measured" stitch backend: tiles placed from their image content, blended at mid-overlap
"""Stage positions alone leave neighbouring CQ1 tiles several pixels apart, by an amount that differs from join to
join, so every join is measured:

1. For each pair of neighbouring tiles, high-pass maximum projections (a few planes) of the two facing edge strips are
   matched by normalised cross-correlation. A join is accepted when the channels agree; otherwise it gets the typical
   offset of the other joins in the same direction, or the stage position if there is none.
2. Tile positions are solved from all joins together (least squares). One position per tile, for every plane and
   channel.
3. Blending: every output pixel is a weighted mean of raw pixels. The weight changes from one tile to the other in the
   middle of the overlap, over `blend_px` pixels. Nothing corrects brightness differences between tiles.

The matching and the position solver are ported from the figure-stitching script `stitch_v3.py` (2026-11-18).
"""
from __future__ import annotations

import json
import time
from concurrent.futures import ThreadPoolExecutor
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from tifffile import imread, imwrite

from stitcher_unified import (_resolve_tile_path, _safe_phys, _tile_chunks, _writer_kwargs, stitch_grid_seamless,
                              tqdm)

HIGHPASS_PX = 15          # radius of the local mean removed before matching (suppresses shading)
MAX_PLANES = 8            # planes in the projection used for matching
SEARCH_ALONG_PX = 60      # largest offset searched along the join (also at most a third of the overlap)
SEARCH_ACROSS_PX = 40     # largest offset searched across it
MIN_WINDOW_PX = 32        # overlaps leaving less than this to compare are not measured
MIN_NCC = 0.5             # several channels: a channel counts from this correlation on ...
AGREE_PX = 3              # ... and the counting channels must agree within this many pixels
SINGLE_NCC = 0.8          # one channel: stricter, and it must also agree with the other joins of its direction
RESIDUAL_WARN_PX = 3
DEFAULT_BLEND_PX = 128
READ_THREADS = 8          # tiles read at the same time (much faster on network drives)

Position = Dict[int, Tuple[int, int]]        # FieldIndex -> (x, y) in mosaic pixels
Pair = Tuple[int, int, str]                  # (field a, field b, "x": b right of a | "y": b below a)


# ------------------------------------------------------------------ matching
def _box_blur(img: np.ndarray, r: int) -> np.ndarray:
    """Mean filter of size (2r+1)^2 via an integral image."""
    p = np.pad(img.astype(np.float64), r + 1, mode="edge")
    c = p.cumsum(0).cumsum(1)
    k = 2 * r + 1
    out = c[k:, k:] - c[:-k, k:] - c[k:, :-k] + c[:-k, :-k]
    return (out[: img.shape[0], : img.shape[1]] / (k * k)).astype(np.float32)


def _highpass(img: np.ndarray) -> np.ndarray:
    img = img.astype(np.float32)
    return img - _box_blur(img, HIGHPASS_PX)


def _fast_size(n: int) -> int:
    """Smallest length >= n with no prime factor above 5 (fast to Fourier transform)."""
    while True:
        m = n
        for p in (2, 3, 5):
            while m % p == 0:
                m //= p
        if m == 1:
            return n
        n += 1


def _box_sums(img: np.ndarray, r0: np.ndarray, r1: np.ndarray, c0: np.ndarray, c1: np.ndarray) -> np.ndarray:
    """Sum of img[r0:r1, c0:c1] for arrays of limits, from a summed-area table."""
    table = np.zeros((img.shape[0] + 1, img.shape[1] + 1))
    table[1:, 1:] = img.cumsum(0).cumsum(1)
    return table[r1, c1] - table[r0, c1] - table[r1, c0] + table[r0, c0]


def _correlation_map(a: np.ndarray, b: np.ndarray, ox: np.ndarray, dy: np.ndarray) -> np.ndarray:
    """Normalised cross-correlation of the overlap for every b origin (ox, dy) relative to a's origin (b to the right
    of a, 0 < ox < width). Exact: the sum of products comes from one Fourier transform, the sums and sums of squares
    over each overlap rectangle from summed-area tables (J. P. Lewis, "Fast Normalized Cross-Correlation",
    https://scribblethink.org/Work/nvisionInterface/nip.html; summed-area tables: F. C. Crow, Computer Graphics
    18(3), 1984, https://doi.org/10.1145/800031.808600). Returns an array of shape (len(ox), len(dy))."""
    h, w = a.shape
    a = a.astype(np.float64)
    b = b.astype(np.float64)
    a -= a.mean()  # keeps the sums small; does not change any correlation
    b -= b.mean()
    size = (_fast_size(h + int(np.abs(dy).max())), _fast_size(w + int(ox.max())))
    spectrum = np.fft.rfft2(a, size) * np.conj(np.fft.rfft2(b, size))
    products = np.fft.irfft2(spectrum, size)[np.ix_(dy % size[0], ox)].T
    x, y = ox[:, None], dy[None, :]
    count = (w - x) * (h - np.abs(y))
    a_rows = (np.maximum(y, 0), h + np.minimum(y, 0))
    b_rows = (np.maximum(-y, 0), h + np.minimum(-y, 0))
    sum_a = _box_sums(a, a_rows[0], a_rows[1], x, np.full_like(x, w))
    sum_b = _box_sums(b, b_rows[0], b_rows[1], np.zeros_like(x), w - x)
    var_a = _box_sums(a * a, a_rows[0], a_rows[1], x, np.full_like(x, w)) - sum_a ** 2 / count
    var_b = _box_sums(b * b, b_rows[0], b_rows[1], np.zeros_like(x), w - x) - sum_b ** 2 / count
    usable = (var_a > 1e-6 * count) & (var_b > 1e-6 * count)  # a flat overlap has no correlation
    out = np.full(products.shape, -1.0)
    out[usable] = ((products - sum_a * sum_b / count)[usable]) / np.sqrt(var_a[usable] * var_b[usable])
    return out


def _measure_pair(a: np.ndarray, b: np.ndarray, nominal_ox: int, nominal_dy: int, r_along: int,
                  r_across: int) -> Tuple[float, int, int]:
    """(ncc, offset along, offset across) of strip b relative to its nominal place: the best of all shifts within
    +-r_along / +-r_across pixels."""
    h, w = a.shape
    ox = np.arange(max(nominal_ox - r_along, 1), min(nominal_ox + r_along, w - 8) + 1)
    dy = np.arange(max(nominal_dy - r_across, 1 - h), min(nominal_dy + r_across, h - 1) + 1)
    if ox.size == 0 or dy.size == 0:
        return -1.0, 0, 0
    scores = _correlation_map(a, b, ox, dy)
    i, j = np.unravel_index(int(np.argmax(scores)), scores.shape)
    return float(scores[i, j]), int(ox[i] - nominal_ox), int(dy[j] - nominal_dy)


# ------------------------------------------------------------------ layout
def stage_positions(g: pd.DataFrame) -> Position:
    """Tile positions in pixels from the stage coordinates, one per field, origin at the top left."""
    if not {"PixelX_stage", "PixelY_stage", "FieldIndex"}.issubset(g.columns):
        raise KeyError("measured placement needs PixelX_stage, PixelY_stage and FieldIndex in the metadata CSV.")
    first = g.drop_duplicates("FieldIndex")
    x0, y0 = first["PixelX_stage"].min(), first["PixelY_stage"].min()
    return {int(r["FieldIndex"]): (int(round(r["PixelX_stage"] - x0)), int(round(r["PixelY_stage"] - y0)))
            for _, r in first.iterrows()}


def neighbour_pairs(stage: Position, tw: int, th: int) -> List[Pair]:
    """Each tile with its nearest right-hand and nearest lower neighbour, where the two overlap or touch."""
    pairs: List[Pair] = []
    for a, (xa, ya) in sorted(stage.items()):
        # "next to" = at least a fifth of a tile further along (overlaps up to 80 %), less than a fifth sideways
        right = [(xb - xa, abs(yb - ya), b) for b, (xb, yb) in stage.items()
                 if tw / 5 < xb - xa <= tw + 8 and abs(yb - ya) < th / 5]
        below = [(yb - ya, abs(xb - xa), b) for b, (xb, yb) in stage.items()
                 if th / 5 < yb - ya <= th + 8 and abs(xb - xa) < tw / 5]
        if right:
            pairs.append((a, min(right)[2], "x"))
        if below:
            pairs.append((a, min(below)[2], "y"))
    return pairs


def measure_joins(g: pd.DataFrame, image_root: str, stage: Position, pairs: List[Pair], tw: int, th: int,
                  loader: Callable = imread, progress_tag: Optional[str] = None) -> Dict[Pair, dict]:
    """For every pair: {"skipped": reason} or {"channels": {c: (ncc, along, across)}}."""
    zs = sorted(g["Z"].unique())
    use = [zs[i] for i in np.unique(np.linspace(0, len(zs) - 1, min(MAX_PLANES, len(zs))).astype(int))]
    result: Dict[Pair, dict] = {}
    plan = {}
    for pair in pairs:
        a, b, axis = pair
        i, size = (0, tw) if axis == "x" else (1, th)
        along = stage[b][i] - stage[a][i]
        across = stage[b][1 - i] - stage[a][1 - i]
        overlap = size - along
        r = min(SEARCH_ALONG_PX, max(overlap // 3, 0))
        if overlap - r < MIN_WINDOW_PX:
            result[pair] = {"skipped": f"overlap of {overlap} px is too small to measure"}
            continue
        strip = min(overlap + r + 8, size)
        plan[pair] = (r, strip, along - (size - strip), across)
        result[pair] = {"channels": {}}

    channels = sorted(g["C"].unique())
    for number, c in enumerate(channels, 1):
        if progress_tag is not None and plan:  # this step takes minutes on large grids
            print(f"{progress_tag}   measuring {len(plan)} joins, channel {number} of {len(channels)} ...", flush=True)
        strips: Dict[Tuple[Pair, int], np.ndarray] = {}
        for field in sorted({f for pair in plan for f in pair[:2]}):
            sub = g[(g["C"] == c) & (g["FieldIndex"] == field) & (g["Z"].isin(use))]
            mip = None
            with ThreadPoolExecutor(READ_THREADS) as pool:
                for img in pool.map(loader, [_resolve_tile_path(name, image_root) for name in sub["TiffFile"]]):
                    mip = img if mip is None else np.maximum(mip, img)
            if mip is None:
                continue
            for pair, (_, s, _, _) in plan.items():  # keep only the facing edge strips (columns = along the join)
                a, b, axis = pair
                view = mip if axis == "x" else mip.T
                if field == a:
                    strips[(pair, 0)] = view[:, -s:].copy()
                if field == b:
                    strips[(pair, 1)] = view[:, :s].copy()
        for pair, (r, _, nominal_ox, nominal_dy) in plan.items():
            if (pair, 0) in strips and (pair, 1) in strips:
                result[pair]["channels"][int(c)] = _measure_pair(
                    _highpass(strips[(pair, 0)]), _highpass(strips[(pair, 1)]),
                    nominal_ox, nominal_dy, r, SEARCH_ACROSS_PX)
    return result


def decide_joins(measured: Dict[Pair, dict], n_channels: int) -> Dict[Pair, dict]:
    """Accept or replace every join. Returns per pair {status, reason, along, across}."""
    out: Dict[Pair, dict] = {}
    candidates: Dict[Pair, Tuple[float, float]] = {}
    for pair, m in measured.items():
        if "skipped" in m:
            out[pair] = {"reason": m["skipped"]}
            continue
        values = list(m["channels"].values())
        if n_channels > 1:
            good = [v for v in values if v[0] >= MIN_NCC]
            if len(good) < 2:
                out[pair] = {"reason": f"fewer than 2 channels reach a correlation of {MIN_NCC}"}
                continue
            along, across = [v[1] for v in good], [v[2] for v in good]
            if max(along) - min(along) > AGREE_PX or max(across) - min(across) > AGREE_PX:
                out[pair] = {"reason": f"channels disagree by more than {AGREE_PX} px"}
                continue
            out[pair] = {"status": "measured", "along": float(np.median(along)), "across": float(np.median(across))}
        elif values and values[0][0] >= SINGLE_NCC:
            candidates[pair] = (values[0][1], values[0][2])
        else:
            out[pair] = {"reason": f"correlation below {SINGLE_NCC} (single channel)"}
    for axis in ("x", "y"):  # single channel: the join must agree with the others of its direction
        same = {p: v for p, v in candidates.items() if p[2] == axis}
        if not same:
            continue
        med = np.median(np.array(list(same.values()), float), axis=0)
        for pair, v in same.items():
            if len(same) >= 2 and abs(v[0] - med[0]) <= AGREE_PX and abs(v[1] - med[1]) <= AGREE_PX:
                out[pair] = {"status": "measured", "along": float(v[0]), "across": float(v[1])}
            else:
                out[pair] = {"reason": "single channel, and no second join of this direction agrees with it"}
    for axis in ("x", "y"):  # replace what could not be measured
        good = [(d["along"], d["across"]) for p, d in out.items() if p[2] == axis and d.get("status") == "measured"]
        med = np.median(np.array(good, float), axis=0) if good else None
        for pair, d in out.items():
            if pair[2] == axis and d.get("status") != "measured":
                if med is None:
                    d.update(status="stage", along=0.0, across=0.0)
                else:
                    d.update(status="typical", along=float(med[0]), across=float(med[1]))
    return out


def solve_positions(stage: Position, constraints: Dict[Tuple[int, int], Tuple[float, float, float]]) -> Position:
    """Least squares for constraints[(a, b)] = (dx, dy, weight), the wanted position of b minus that of a.
    The first field is anchored at its stage position; a weak stage prior holds tiles without any constraint."""
    fields = sorted(stage)
    index = {f: i for i, f in enumerate(fields)}
    solved = []
    for axis in (0, 1):
        rows, rhs = [], []
        for (a, b), value in constraints.items():
            v = np.zeros(len(fields))
            v[index[b]], v[index[a]] = 1, -1
            weight = np.sqrt(value[2])
            rows.append(v * weight)
            rhs.append(value[axis] * weight)
        v = np.zeros(len(fields))
        v[0] = 1
        rows.append(v)
        rhs.append(stage[fields[0]][axis])
        for f in fields:
            v = np.zeros(len(fields))
            v[index[f]] = 1e-3
            rows.append(v)
            rhs.append(1e-3 * stage[f][axis])
        solved.append(np.linalg.lstsq(np.array(rows), np.array(rhs), rcond=None)[0])
    x0, y0 = solved[0].min(), solved[1].min()
    return {f: (int(round(solved[0][index[f]] - x0)), int(round(solved[1][index[f]] - y0))) for f in fields}


def place_tiles(g: pd.DataFrame, image_root: str, tw: int, th: int, tag: str = "",
                verbose: bool = True) -> Tuple[Position, dict]:
    """Measured tile positions for one (well, grid) and a report of how each join was placed."""
    t0 = time.perf_counter()
    stage = stage_positions(g)
    pairs = neighbour_pairs(stage, tw, th)
    n_channels = int(g["C"].nunique())
    measured = measure_joins(g, image_root, stage, pairs, tw, th, progress_tag=tag if verbose else None)
    decided = decide_joins(measured, n_channels)
    constraints = {}
    for pair, d in decided.items():
        a, b, axis = pair
        base = (stage[b][0] - stage[a][0], stage[b][1] - stage[a][1])
        dx, dy = (d["along"], d["across"]) if axis == "x" else (d["across"], d["along"])
        constraints[(a, b)] = (base[0] + dx, base[1] + dy, 1.0 if d["status"] == "measured" else 0.1)
    positions = solve_positions(stage, constraints)

    joins, warnings = [], []
    for pair, d in sorted(decided.items()):
        a, b, axis = pair
        want = constraints[(a, b)]
        residual = float(np.hypot(positions[b][0] - positions[a][0] - want[0],
                                  positions[b][1] - positions[a][1] - want[1]))
        entry = dict(fields=[a, b], direction="left-right" if axis == "x" else "top-bottom", status=d["status"],
                     offset_along_px=d["along"], offset_across_px=d["across"], reason=d.get("reason", ""),
                     channels={str(c + 1): [round(v[0], 3), v[1], v[2]]
                               for c, v in measured[pair].get("channels", {}).items()},
                     residual_px=round(residual, 1))
        joins.append(entry)
        if d["status"] == "measured" and residual > RESIDUAL_WARN_PX:
            warnings.append(f"fields {a}-{b}: measured offset and solved positions differ by {residual:.1f} px")
    counts = {s: sum(1 for j in joins if j["status"] == s) for s in ("measured", "typical", "stage")}
    if verbose:
        print(f"{tag} tile placement: {len(joins)} joins - {counts['measured']} measured, "
              f"{counts['typical']} given the typical offset, {counts['stage']} left at the stage position "
              f"({time.perf_counter() - t0:.0f} s)")
        for j in joins:
            if j["status"] != "measured":
                print(f"{tag}   fields {j['fields'][0]}-{j['fields'][1]} ({j['direction']}): {j['status']} - "
                      f"{j['reason']}")
        for w in warnings:
            print(f"⚠️ {tag} {w}")
    report = dict(method="measured", joins=joins, counts=counts, warnings=warnings,
                  stage_positions_px={str(f): list(p) for f, p in stage.items()})
    return positions, report


# ------------------------------------------------------------------ blending
def blend_width(positions: Position, pairs: List[Pair], tw: int, th: int, requested: int) -> int:
    """`requested`, but not more than the smallest overlap between neighbouring tiles (a wider blend would leave a
    step at the border of that overlap)."""
    overlaps = [(tw - (positions[b][0] - positions[a][0])) if axis == "x" else (th - (positions[b][1] - positions[a][1]))
                for a, b, axis in pairs]
    positive = [o for o in overlaps if o > 0]
    return max(1, min([int(requested)] + positive))


def blend_plane(paths: Dict[int, str], positions: Position, tw: int, th: int, width: int, height: int,
                blend_px: int, dtype, loader: Callable = imread, band: int = 1024) -> np.ndarray:
    """One mosaic plane. Weight of a tile at a pixel: clip(1 - (deepest - own depth) / blend_px, 0, 1), where depth is
    the distance to the tile's own edge and `deepest` the largest depth among the tiles covering the pixel. Works in
    horizontal bands, so the only full-mosaic array is the output plane."""
    plane = np.zeros((height, width), dtype=dtype)
    ax = np.minimum(np.arange(tw), np.arange(tw)[::-1]).astype(np.float32)
    ay = np.minimum(np.arange(th), np.arange(th)[::-1]).astype(np.float32)
    cache: Dict[int, np.ndarray] = {}
    for y0 in range(0, height, band):
        y1 = min(y0 + band, height)
        for f in [f for f in cache if positions[f][1] + th <= y0]:
            del cache[f]
        here = []
        for f in sorted(paths):
            x, y = positions[f]
            if y < y1 and y + th > y0:
                r0, r1 = max(y0, y), min(y1, y + th)
                here.append((f, x, r0 - y, r1 - y, slice(r0 - y0, r1 - y0)))
        deepest = np.full((y1 - y0, width), -1.0, np.float32)
        for f, x, t0, t1, rows in here:
            view = deepest[rows, x:x + tw]
            np.maximum(view, np.minimum(ax[None, :], ay[t0:t1, None]), out=view)
        total = np.zeros((y1 - y0, width), np.float32)
        weight = np.zeros((y1 - y0, width), np.float32)
        wanted = [f for f, *_ in here if f not in cache]
        if wanted:
            with ThreadPoolExecutor(READ_THREADS) as pool:
                cache.update(zip(wanted, pool.map(loader, [paths[f] for f in wanted])))
        for f, x, t0, t1, rows in here:
            w =np.clip(1 - (deepest[rows, x:x + tw] - np.minimum(ax[None, :], ay[t0:t1, None])) / blend_px, 0, 1)
            total[rows, x:x + tw] += w * cache[f][t0:t1]
            weight[rows, x:x + tw] += w
        total /= np.maximum(weight, 1e-6)
        if np.issubdtype(plane.dtype, np.integer):
            np.rint(total, out=total)
            np.clip(total, 0, np.iinfo(plane.dtype).max, out=total)
        plane[y0:y1] = total
    return plane


# ------------------------------------------------------------------ backend
def _write_report(output_path: str, report: dict) -> None:
    """Write <name>_placement.json next to the stitched file."""
    stem = output_path[:-len(".ome.tif")] if output_path.endswith(".ome.tif") else output_path
    with open(stem + "_placement.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=1)


def stitch_grid_measured(
    df: pd.DataFrame, image_root: str, output_path: str, well_id: str, grid_index: int, *,
    placement_df: Optional[pd.DataFrame] = None, positions: Optional[Position] = None,
    blend_px: int = DEFAULT_BLEND_PX, verbose: bool = True,
    compression: Optional[str] = "zlib", compression_level: int = 6, predictor: bool = True,
    tile: Optional[Tuple[int, int]] = (512, 512), bigtiff: bool = True, pyramids: bool = False,
    # stage method, used only when a grid falls back to it
    default_overlap_fraction: float = 0.01, min_feather_px: int = 4, max_feather_px: int = 16,
    max_snap_px: int = 6, lattice_snap_px: int = 0, snap_axis: str = "both", refine_alignment: bool = False,
    refine_max_shift_px: int = 6, refine_axis: str = "both", overlap_gain_match: bool = False,
    gain_clip_low: float = 0.92, gain_clip_high: float = 1.08, gain_match_channels: Optional[list] = None,
) -> None:
    """Stitch one (well, grid). `placement_df`: metadata used to measure positions (all planes, even when only one is
    stitched). `positions`: given tile positions instead of measured ones (for tests).

    If any join can neither be measured nor take the typical offset of its direction (it would stay at the stage
    position, several pixels off), the whole grid is stitched with the stage method instead: blending misaligned
    tiles over a wide band doubles structures."""
    t0 = time.perf_counter()
    g = df[(df["WellID"] == well_id) & (df["GridIndex"] == grid_index)].copy()
    if g.empty:
        print(f"No data for Well {well_id}, Grid {grid_index}. Skipping.")
        return
    img0 = imread(_resolve_tile_path(g["TiffFile"].iloc[0], image_root))
    th, tw = img0.shape
    dtype = img0.dtype
    psx = _safe_phys(g["PxSizeX_um"].iloc[0])
    psy = _safe_phys(g["PxSizeY_um"].iloc[0])
    psz = _safe_phys(g.get("PxSizeZ_um", pd.Series([1.0])).iloc[0])
    tag = f"[{well_id} G{grid_index}]"

    source = g if placement_df is None else placement_df[(placement_df["WellID"] == well_id)
                                                         & (placement_df["GridIndex"] == grid_index)]
    if positions is None:
        if verbose:
            print(f"{tag} measured placement (the feather, snap and refine settings apply to the stage method only)")
        positions, report = place_tiles(source, image_root, tw, th, tag, verbose)
        if report["counts"]["stage"]:
            reason = (f"{report['counts']['stage']} of {len(report['joins'])} joins could not be placed from the "
                      f"images")
            print(f"⚠️ {tag} {reason}, so this grid is stitched with the stage method (positions from the "
                  f"microscope, narrow blend).")
            stitch_grid_seamless(
                df, image_root, output_path, well_id, grid_index, placement="stage", verbose=verbose,
                compression=compression, compression_level=compression_level, predictor=predictor, tile=tile,
                bigtiff=bigtiff, pyramids=pyramids, default_overlap_fraction=default_overlap_fraction,
                min_feather_px=min_feather_px, max_feather_px=max_feather_px, max_snap_px=max_snap_px,
                lattice_snap_px=lattice_snap_px, snap_axis=snap_axis, refine_alignment=refine_alignment,
                refine_max_shift_px=refine_max_shift_px, refine_axis=refine_axis,
                overlap_gain_match=overlap_gain_match, gain_clip_low=gain_clip_low, gain_clip_high=gain_clip_high,
                gain_match_channels=gain_match_channels)
            report.update(method="stage (fallback)", fallback_reason=reason)
            _write_report(output_path, report)
            return
    else:
        report = dict(method="given", joins=[], counts={}, warnings=[])
    pairs = neighbour_pairs(stage_positions(source), tw, th)
    width_px = blend_width(positions, pairs, tw, th, blend_px)
    max_x = max(x for x, _ in positions.values()) + tw
    max_y = max(y for _, y in positions.values()) + th
    if verbose:
        print(f"{tag} mosaic size: {max_x}x{max_y} px; tile {tw}x{th}; blend {width_px} px"
              + (f" (limited by the smallest overlap, asked for {blend_px})" if width_px < blend_px else ""))

    z_slices = sorted(g["Z"].unique())
    channels = sorted(g["C"].unique())
    out_shape = (len(z_slices), len(channels), max_y, max_x)
    progress = tqdm(total=len(z_slices) * len(channels), desc=f"Blending planes {well_id} G{grid_index}", unit="plane")

    def _planes():
        for z in z_slices:
            for c in channels:
                rows = g[(g["Z"] == z) & (g["C"] == c)]
                paths = {int(r["FieldIndex"]): _resolve_tile_path(r["TiffFile"], image_root)
                         for _, r in rows.iterrows()}
                plane = blend_plane(paths, positions, tw, th, max_x, max_y, width_px, dtype)
                progress.update(1)
                yield plane

    # one name per channel index (two channels may carry the same name)
    ch_names = [g.loc[g["C"] == c, "ChannelName"].iloc[0] for c in channels]
    metadata = {"axes": "ZCYX", "PhysicalSizeX": psx, "PhysicalSizeY": psy, "PhysicalSizeZ": psz,
                "Channel": {"Name": ch_names}}
    try:
        imwrite(output_path, _tile_chunks(_planes(), tile, max_y, max_x), shape=out_shape, dtype=dtype,
                metadata=metadata,
                **_writer_kwargs(compression, compression_level, predictor, tile, bigtiff, pyramids))
    finally:
        progress.close()
    report.update(positions_px={str(f): list(p) for f, p in positions.items()}, tile_px=[tw, th],
                  mosaic_px=[max_x, max_y], blend_px=width_px)
    _write_report(output_path, report)
    print(f" Saved OME-TIFF: {output_path}  Shape: {out_shape} (Z,C,Y,X)  total={time.perf_counter() - t0:.2f}s")
