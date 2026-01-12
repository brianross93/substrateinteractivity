import argparse
import os
from pathlib import Path

import numpy as np


def _count_local_minima(y: np.ndarray, prominence: float) -> int:
    """
    Very simple dip counter: counts points that are lower than both neighbors by at least `prominence`.
    """
    if y.size < 3:
        return 0
    a = y[:-2]
    b = y[1:-1]
    c = y[2:]
    dips = (b < a) & (b < c) & ((np.minimum(a, c) - b) >= prominence)
    return int(np.sum(dips))


def _top_dips(wl: np.ndarray, y: np.ndarray, prominence: float, top_n: int) -> list[tuple[float, float, float]]:
    """
    Returns list of (wavelength, T_at_dip, dip_depth) sorted by depth desc.
    dip_depth = min(neighbors) - dip.
    """
    if y.size < 3 or top_n <= 0:
        return []
    a = y[:-2]
    b = y[1:-1]
    c = y[2:]
    depth = np.minimum(a, c) - b
    dips = (b < a) & (b < c) & (depth >= prominence)
    idxs = np.where(dips)[0] + 1
    items = [(float(wl[i]), float(y[i]), float(depth[i - 1])) for i in idxs]
    items.sort(key=lambda t: t[2], reverse=True)
    return items[:top_n]


def summarize_curve(name: str, wl: np.ndarray, t: np.ndarray) -> dict:
    idx_min = int(np.argmin(t))
    idx_max = int(np.argmax(t))
    return {
        "name": name,
        "wl_minT": float(wl[idx_min]),
        "minT": float(t[idx_min]),
        "wl_maxT": float(wl[idx_max]),
        "maxT": float(t[idx_max]),
        "meanT": float(np.mean(t)),
        "medianT": float(np.median(t)),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", default="outputs/cantor_transmission_wideband.npz", help="Combined NPZ file from meep_cantor.py")
    ap.add_argument("--dip-prominence", type=float, default=0.02, help="Minimum dip depth to count as a dip")
    ap.add_argument("--top-dips", type=int, default=5, help="Print the top N deepest dips per curve")
    args = ap.parse_args()

    npz_path = Path(args.npz)
    if not npz_path.exists():
        raise SystemExit(f"Missing: {npz_path}")

    data = np.load(npz_path)
    wl = data["wavelength_um"].astype(np.float64)
    levels = data["levels"].astype(np.int32).tolist()
    transmissions = data["transmissions"].astype(np.float64)  # (L, N)

    print(f"NPZ: {os.path.abspath(npz_path)}")
    print(f"Wavelength range: {wl.min():.3f} .. {wl.max():.3f} um (N={wl.size})")
    print(f"Levels in file: {levels}")
    print(f"Dip counter prominence: {args.dip_prominence:.3f}")
    print("")

    for i, L in enumerate(levels):
        t = transmissions[i]
        s = summarize_curve(f"Cantor L={L}", wl, t)
        dips = _count_local_minima(t, prominence=float(args.dip_prominence))
        print(
            f"{s['name']}: "
            f"minT={s['minT']:.3f} @ {s['wl_minT']:.3f}um, "
            f"maxT={s['maxT']:.3f} @ {s['wl_maxT']:.3f}um, "
            f"meanT={s['meanT']:.3f}, medianT={s['medianT']:.3f}, "
            f"dips~{dips}"
        )
        top = _top_dips(wl, t, prominence=float(args.dip_prominence), top_n=int(args.top_dips))
        if top:
            for j, (wlj, tj, dj) in enumerate(top, start=1):
                print(f"  dip#{j}: wl={wlj:.3f}um  T={tj:.3f}  depth={dj:.3f}")

    # Optional slab baseline (if present): accept any *_slab_match_L*.npz
    slab = None
    for cand in sorted(npz_path.parent.glob(npz_path.stem + "_slab_match_L*.npz")):
        slab = cand
        break

    if slab is not None and slab.exists():
        sdata = np.load(slab)
        swl = sdata["wavelength_um"].astype(np.float64)
        st = sdata["transmission"].astype(np.float64)
        s = summarize_curve("Slab baseline", swl, st)
        dips = _count_local_minima(st, prominence=float(args.dip_prominence))
        print("")
        print(
            f"{s['name']}: "
            f"minT={s['minT']:.3f} @ {s['wl_minT']:.3f}um, "
            f"maxT={s['maxT']:.3f} @ {s['wl_maxT']:.3f}um, "
            f"meanT={s['meanT']:.3f}, medianT={s['medianT']:.3f}, "
            f"dips~{dips}"
        )
        top = _top_dips(swl, st, prominence=float(args.dip_prominence), top_n=int(args.top_dips))
        if top:
            for j, (wlj, tj, dj) in enumerate(top, start=1):
                print(f"  dip#{j}: wl={wlj:.3f}um  T={tj:.3f}  depth={dj:.3f}")
    else:
        print("")
        print("Slab baseline NPZ not found (expected something like *_slab_match_L*.npz).")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

