import argparse
import csv
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def _img_page(pdf: PdfPages, title: str, img_path: Path) -> None:
    fig = plt.figure(figsize=(8.5, 11), constrained_layout=True)
    ax = fig.add_subplot(111)
    ax.set_title(title, fontsize=14)
    ax.axis("off")
    if img_path.exists():
        img = plt.imread(str(img_path))
        ax.imshow(img)
    else:
        ax.text(0.5, 0.5, f"Missing: {img_path.name}", ha="center", va="center", fontsize=12)
    pdf.savefig(fig)
    plt.close(fig)


def _table_page(pdf: PdfPages, title: str, csv_path: Path, max_rows: int = 40) -> None:
    fig = plt.figure(figsize=(8.5, 11), constrained_layout=True)
    ax = fig.add_subplot(111)
    ax.set_title(title, fontsize=14)
    ax.axis("off")

    if not csv_path.exists():
        ax.text(0.5, 0.5, f"Missing: {csv_path.name}", ha="center", va="center", fontsize=12)
        pdf.savefig(fig)
        plt.close(fig)
        return

    rows: list[list[str]] = []
    with open(csv_path, "r", encoding="utf-8") as f:
        r = csv.reader(f)
        for i, row in enumerate(r):
            if i >= max_rows:
                rows.append(["… truncated …"])
                break
            rows.append(row)

    # Matplotlib requires a rectangular table; pad short rows.
    cols = max((len(r) for r in rows), default=1)
    padded = [r + [""] * (cols - len(r)) for r in rows]
    table = ax.table(cellText=padded, loc="upper left", cellLoc="left")
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.0, 1.2)
    pdf.savefig(fig)
    plt.close(fig)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", default="outputs", help="Directory with generated PNG/CSV outputs.")
    p.add_argument("--outfile", default=None, help="Output PDF path (defaults to <outdir>/phase1_report.pdf).")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    outfile = Path(args.outfile) if args.outfile else (outdir / "phase1_report.pdf")

    with PdfPages(str(outfile)) as pdf:
        # Cover / summary
        fig = plt.figure(figsize=(8.5, 11), constrained_layout=True)
        ax = fig.add_subplot(111)
        ax.axis("off")
        ax.set_title("Phase 1 Simulation Report", fontsize=18, pad=20)
        ax.text(
            0.02,
            0.95,
            "This PDF bundles the latest plots + raw CSV tables produced by the headless runs.\n"
            "Regenerate everything via README commands, then re-run this script.\n",
            va="top",
            fontsize=11,
        )
        pdf.savefig(fig)
        plt.close(fig)

        # Plots
        _img_page(pdf, "Kuramoto (r(t) and final phases)", outdir / "kuramoto.png")
        _img_page(pdf, "Kuramoto animation (GIF) – not embedded; see outputs/kuramoto.gif", outdir / "kuramoto.png")
        _img_page(pdf, "Rectangular tunneling (|psi|^2 and V)", outdir / "tunneling.png")
        _img_page(pdf, "CHSH (single run)", outdir / "chsh.png")
        _img_page(pdf, "CHSH noise sweep", outdir / "chsh_noise_sweep.png")
        _img_page(pdf, "Retro (order dependence + time-series)", outdir / "retro.png")
        _img_page(pdf, "Retro spectrum", outdir / "retro_spectrum.png")
        _img_page(pdf, "Quantum switch (CPTP) interference / visibility", outdir / "retro_switch.png")
        _img_page(pdf, "Retro visibility vs theta", outdir / "retro_vis_vs_theta.png")
        _img_page(pdf, "Retro visibility vs phi", outdir / "retro_vis_vs_phi.png")
        _img_page(pdf, "Cantor tunneling sweep + fit", outdir / "cantor_t_sweep.png")

        # CSV tables
        _table_page(pdf, "CSV: CHSH noise sweep", outdir / "chsh_noise_sweep.csv")
        _table_page(pdf, "CSV: Retro visibility vs theta", outdir / "retro_vis_vs_theta.csv")
        _table_page(pdf, "CSV: Cantor tunneling sweep + fit", outdir / "cantor_t_sweep.csv")

    print(f"Wrote {os.path.abspath(outfile)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

