import argparse
import csv
from pathlib import Path

import numpy as np


def load_two_col_csv(path: Path) -> tuple[np.ndarray, np.ndarray]:
    wl = []
    val = []
    with open(path, "r", encoding="utf-8") as f:
        r = csv.reader(f)
        for row in r:
            if not row or row[0].strip().startswith("#"):
                continue
            if row[0].lower().strip() in ("wavelength", "lambda", "wl"):
                continue
            if len(row) < 2:
                continue
            wl.append(float(row[0]))
            val.append(float(row[1]))
    return np.array(wl, dtype=np.float64), np.array(val, dtype=np.float64)


def top_dips(wl: np.ndarray, t: np.ndarray, prominence: float, top_n: int) -> list[tuple[float, float, float]]:
    if t.size < 3 or top_n <= 0:
        return []
    a = t[:-2]
    b = t[1:-1]
    c = t[2:]
    depth = np.minimum(a, c) - b
    dips = (b < a) & (b < c) & (depth >= prominence)
    idxs = np.where(dips)[0] + 1
    items = [(float(wl[i]), float(t[i]), float(depth[i - 1])) for i in idxs]
    items.sort(key=lambda x: x[2], reverse=True)
    return items[:top_n]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True, help="CSV with columns: wavelength_um,power")
    ap.add_argument("--sample", required=True, help="CSV with columns: wavelength_um,power")
    ap.add_argument("--slab", required=False, help="Optional slab-control CSV")
    ap.add_argument("--prominence", type=float, default=0.02)
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--out", default="outputs/measurement_normalized.csv")
    args = ap.parse_args()

    wl_r, p_r = load_two_col_csv(Path(args.ref))
    wl_s, p_s = load_two_col_csv(Path(args.sample))
    if wl_r.size != wl_s.size or np.max(np.abs(wl_r - wl_s)) > 1e-9:
        raise SystemExit("ref/sample wavelength grids must match.")

    t_sample = p_s / np.maximum(p_r, 1e-30)

    t_slab = None
    if args.slab:
        wl_b, p_b = load_two_col_csv(Path(args.slab))
        if wl_b.size != wl_r.size or np.max(np.abs(wl_b - wl_r)) > 1e-9:
            raise SystemExit("slab wavelength grid must match.")
        t_slab = p_b / np.maximum(p_r, 1e-30)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        if t_slab is None:
            w.writerow(["wavelength_um", "T_sample"])
            for wl, t in zip(wl_r.tolist(), t_sample.tolist(), strict=False):
                w.writerow([wl, t])
        else:
            w.writerow(["wavelength_um", "T_sample", "T_slab"])
            for wl, ts, tb in zip(wl_r.tolist(), t_sample.tolist(), t_slab.tolist(), strict=False):
                w.writerow([wl, ts, tb])

    print(f"Wrote {out_path}")
    print("")
    print("Top dips (sample):")
    for i, (wlv, tv, dv) in enumerate(top_dips(wl_r, t_sample, prominence=float(args.prominence), top_n=int(args.top)), start=1):
        print(f"  dip#{i}: wl={wlv:.3f}um  T={tv:.3f}  depth={dv:.3f}")

    if t_slab is not None:
        print("")
        print("Top dips (slab):")
        for i, (wlv, tv, dv) in enumerate(top_dips(wl_r, t_slab, prominence=float(args.prominence), top_n=int(args.top)), start=1):
            print(f"  dip#{i}: wl={wlv:.3f}um  T={tv:.3f}  depth={dv:.3f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

