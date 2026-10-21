"""Shared style for all manuscript figures. Run scripts from /shared/rc/malariarna/MIF_revision_2026."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker  # noqa: F401
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 8, "axes.linewidth": 0.7, "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "ps.fonttype": 42})
OUT = "figures_final"
os.makedirs(OUT, exist_ok=True)
RED, BLUE, GREY = "#C0392B", "#2E6DB4", "#BBBBBB"

def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.eps", format="eps", bbox_inches="tight")
    plt.close(fig)
    print(f"saved {OUT}/{name}.png and .eps")

def tag(ax, s):
    ax.text(-0.2, 1.08, s, transform=ax.transAxes, fontsize=11, fontweight="bold")

def need(path):
    if not os.path.exists(path):
        raise SystemExit(f"Missing input: {path}  (run the analysis script that creates it first)")
    return path
