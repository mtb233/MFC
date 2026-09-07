import glob
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d


def load_bubble_evolution(filename, length_scale, time_scale):
    data = pd.read_csv(filename, sep=r"\s+", skiprows=[1])
    time = data.currentTime * time_scale
    radius = data.radius * length_scale
    return time, radius

def rmse_percent(x_ref, y_ref, x, y, kind="linear"):
    interpolator = interp1d(x, y, kind=kind, fill_value="extrapolate")
    y_interp = interpolator(x_ref)
    rmse = np.sqrt(np.mean((y_ref - y_interp) ** 2))
    return 100 * rmse / np.mean(np.abs(y_ref))


if __name__ == "__main__":
    frmt = "png"

    x0 = 1e-3       # m
    c0 = 1475       # m/s
    t0 = x0/c0      # s

    R0 = 50e-6      # m
    freq = 150e3    # Hz

    # =============================================================================
    # Load testing and golden MFC solutions
    # =============================================================================

    for path in Path("projects/kernel_study/cases/").iterdir():
        if path.is_file():
            continue
        
        files = glob.glob(str(path / "D/lag_bubble_evol_*dat"))
        if files:
            case_name = path.name
            largest_file = max(files, key=os.path.getsize)
            time, radius = load_bubble_evolution(largest_file, x0/R0, t0*1e6)

            time_golden, radius_golden = load_bubble_evolution("projects/kernel_study/post_process/D_golden/lag_bubble_evol_7.dat", x0/R0, t0*1e6)

            # =============================================================================
            # Load analytical solution
            # =============================================================================

            analytical = (pd.read_csv("projects/kernel_study/post_process/analytical_solution.csv", header=None).sort_values(by=0))
            time_offset = 3.24
            time_analytical = analytical[0] * 1e6 / freq + time_offset
            radius_analytical = analytical[1]

            # =============================================================================
            # Error metrics
            # =============================================================================

            golden_error = rmse_percent(time_golden, radius_golden, time, radius)
            analytical_error = rmse_percent(time_analytical, radius_analytical, time, radius)

            # =============================================================================
            # Plot results
            # =============================================================================

            fig, ax = plt.subplots(figsize=(5.5, 4.0))

            ax.plot(time, radius, color="blue", linewidth=1, label="Testing")

            ax.plot(time_golden[::10], radius_golden[::10], "o", markersize=2.5, 
                markerfacecolor="none", markeredgecolor="black", label="Golden")

            ax.plot(time_analytical, radius_analytical, "--", color="black",
                linewidth=1, label="Analytical solution")

            ax.text(15, 1.00, f"RMSE (golden): {golden_error:.2f}%", fontsize='small')
            ax.text(15, 0.90, f"RMSE (analytical): {analytical_error:.2f}%", fontsize='small')

            ax.set_xlim(time_analytical.min(), time_analytical.max())
            ax.set_ylim(0.5, 2.0)

            ax.set_xlabel("Time (μs)")
            ax.set_ylabel(r"$R/R_0$")

            ax.legend(fontsize="small", ncol=3, handlelength=1)

            fig.tight_layout()
            fig.savefig(f"projects/kernel_study/results/{case_name}.{frmt}", dpi=300)
