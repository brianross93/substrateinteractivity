import math
import os
from dataclasses import dataclass

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


@dataclass(frozen=True)
class CantorL3Spec:
    levels: int = 3
    a_um: float = 4.0  # total Cantor window length along propagation
    substrate: str = "Fused silica (SiO2) slide, 25x25x1 mm (typical)"
    concept: str = (
        "1D Cantor multilayer along propagation: alternating high-index material and void/low-index gaps\n"
        "with the Cantor *complement* occupied by the high-index material.\n"
        "Measured in free-space normal incidence through the patterned region."
    )
    target_band: str = "0.6–1.6 µm"
    index_hi: float = 1.5
    index_lo: float = 1.0
    min_feature_nm: float = 4.0e3 / (3.0**3)  # a/3^L in nm, a=4µm
    tolerance_nm: int = 50
    notes: str = (
        "- Start with levels=3 (min feature ~148 nm for a=4 µm).\n"
        "- Baseline control: uniform slab with same total high-index length.\n"
        "- Pattern should extend across the illuminated aperture (e.g., 0.5–2 mm in y).\n"
    )


def cantor_complement_intervals(a_um: float, levels: int) -> list[tuple[float, float]]:
    """
    Returns intervals (start,end) within [-a/2, a/2] that belong to the Cantor complement at given level.
    Constructed analytically via recursive removal.
    """
    intervals = [(-a_um / 2.0, a_um / 2.0)]
    removed: list[tuple[float, float]] = []
    for _ in range(levels):
        new_intervals = []
        for (s, e) in intervals:
            L = e - s
            a1 = s
            b1 = s + L / 3.0
            a2 = s + 2.0 * L / 3.0
            b2 = e
            removed.append((b1, a2))  # middle third removed => complement material
            new_intervals.append((a1, b1))
            new_intervals.append((a2, b2))
        intervals = new_intervals
    # removed may overlap across recursion; merge + sort
    removed.sort()
    merged: list[tuple[float, float]] = []
    for s, e in removed:
        if not merged or s > merged[-1][1]:
            merged.append((s, e))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
    return merged


def main() -> int:
    spec = CantorL3Spec()
    out_pdf = os.path.abspath(os.path.join("outputs", "fab_spec_cantor_L3.pdf"))
    os.makedirs(os.path.dirname(out_pdf), exist_ok=True)

    intervals = cantor_complement_intervals(spec.a_um, spec.levels)
    total_hi = sum(e - s for s, e in intervals)

    with PdfPages(out_pdf) as pdf:
        fig = plt.figure(figsize=(8.5, 11), constrained_layout=True)
        gs = fig.add_gridspec(12, 1)

        ax0 = fig.add_subplot(gs[0:2, 0])
        ax0.axis("off")
        ax0.set_title("Fabrication Spec (MVP): 1D Optical Cantor Barrier (Level 3)", fontsize=16, pad=10)
        ax0.text(
            0.0,
            0.8,
            f"Target band: {spec.target_band}\n"
            f"Model indices: n_hi={spec.index_hi}, n_lo={spec.index_lo}\n"
            f"Cantor window length a = {spec.a_um:.3f} µm, levels={spec.levels}\n"
            f"Min feature ≈ a/3^L = {spec.a_um/(3**spec.levels):.6f} µm ≈ {spec.min_feature_nm:.0f} nm\n"
            f"Total high-index length (within a): {total_hi:.6f} µm\n"
            f"Suggested tolerance: ±{spec.tolerance_nm} nm (feature widths)\n",
            va="top",
            fontsize=11,
        )

        ax1 = fig.add_subplot(gs[2:5, 0])
        ax1.set_title("1D layout along propagation (x) — high-index segments (Cantor complement)", fontsize=12)
        ax1.set_xlim(-spec.a_um / 2.0, spec.a_um / 2.0)
        ax1.set_ylim(0, 1)
        ax1.set_yticks([])
        for s, e in intervals:
            ax1.fill_between([s, e], 0.2, 0.8, color="#2a6fdb", alpha=0.85)
        ax1.axvline(0, color="black", lw=1, alpha=0.3)
        ax1.set_xlabel("x (µm)")
        ax1.text(
            0.02,
            0.05,
            "Blue = high-index material regions; gaps = low-index/void.\n"
            "This is the target geometry used by the Meep 0.6–1.6 µm L=3 run.",
            transform=ax1.transAxes,
            fontsize=9,
            va="bottom",
        )

        ax2 = fig.add_subplot(gs[5:10, 0])
        ax2.axis("off")
        ax2.set_title("Segment table (x-start, x-end, length) [µm]", fontsize=12, pad=6)
        rows = [["#", "x_start", "x_end", "len"]]
        for i, (s, e) in enumerate(intervals, start=1):
            rows.append([str(i), f"{s:.6f}", f"{e:.6f}", f"{(e-s):.6f}"])
        table = ax2.table(cellText=rows, loc="upper left", cellLoc="left")
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1.0, 1.25)

        ax3 = fig.add_subplot(gs[10:12, 0])
        ax3.axis("off")
        ax3.set_title("Fabrication notes", fontsize=12, pad=6)
        ax3.text(0.0, 0.9, f"Substrate: {spec.substrate}\n", fontsize=10, va="top")
        ax3.text(0.0, 0.72, spec.concept, fontsize=10, va="top")
        ax3.text(0.0, 0.25, spec.notes, fontsize=10, va="top")

        pdf.savefig(fig)
        plt.close(fig)

    print(f"Wrote {out_pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

