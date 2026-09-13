import csv
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm

EPS0 = 8.8541878128e-12
K = 1.0 / (4.0 * np.pi * EPS0)
OUTPUT_DIR = Path("outputs")
OUTPUT_DIR.mkdir(exist_ok=True)


def read_float(prompt, default=None):
    while True:
        suffix = f" [{default}]" if default is not None else ""
        value = input(f"{prompt}{suffix}: ").strip()
        if not value and default is not None:
            return float(default)
        try:
            return float(value)
        except ValueError:
            print("Enter a valid number.")


def read_int(prompt, default=None, minimum=None):
    while True:
        suffix = f" [{default}]" if default is not None else ""
        value = input(f"{prompt}{suffix}: ").strip()
        if not value and default is not None:
            return int(default)
        try:
            value = int(value)
            if minimum is not None and int(value) < minimum:
                raise ValueError
            return int(value)
        except ValueError:
            print("Enter a valid integer.")


def charge_presets():
    return {
        "1": [
            (-0.5, 0.0, 1e-9),
            (0.5, 0.0, -1e-9),
        ],
        "2": [
            (-0.5, -0.5, 1e-9),
            (-0.5, 0.5, -1e-9),
            (0.5, -0.5, -1e-9),
            (0.5, 0.5, 1e-9),
        ],
        "3": [
            (-0.6, 0.0, 1e-9),
            (0.0, 0.0, -2e-9),
            (0.6, 0.0, 1e-9),
        ],
        "4": [
            (-0.7, -0.4, 1e-9),
            (0.7, -0.4, -1e-9),
            (-0.7, 0.4, -1e-9),
            (0.7, 0.4, 1e-9),
            (0.0, 0.0, 0.5e-9),
        ],
    }


def read_charges():
    presets = charge_presets()

    print("\nCharge configuration")
    print("1  Electric dipole")
    print("2  Quadrupole")
    print("3  Three-charge system")
    print("4  Five-charge system")
    print("5  Custom configuration")

    choice = input("Select configuration: ").strip()

    if choice in presets:
        return presets[choice]

    count = read_int("Number of point charges", 2, 1)
    charges = []

    for i in range(count):
        print(f"\nCharge {i + 1}")
        x = read_float("x position (m)")
        y = read_float("y position (m)")
        q = read_float("charge (C)")
        charges.append((x, y, q))

    return charges


def compute_field_potential(X, Y, charges):
    ex = np.zeros_like(X, dtype=float)
    ey = np.zeros_like(Y, dtype=float)
    potential = np.zeros_like(X, dtype=float)

    for xq, yq, q in charges:
        dx = X - xq
        dy = Y - yq
        r2 = dx * dx + dy * dy
        r2_safe = np.maximum(r2, 1e-14)
        r = np.sqrt(r2_safe)

        ex += K * q * dx / (r2_safe * r)
        ey += K * q * dy / (r2_safe * r)
        potential += K * q / r

    return ex, ey, potential


def numerical_validation(ex, ey, potential, x, y):
    dV_dy, dV_dx = np.gradient(potential, y, x)
    ex_gradient = -dV_dx
    ey_gradient = -dV_dy

    difference = np.hypot(ex - ex_gradient, ey - ey_gradient)
    field_scale = np.maximum(np.hypot(ex, ey), 1e-14)
    error = difference / field_scale

    return ex_gradient, ey_gradient, error


def orthogonality_analysis(ex, ey, potential, x, y):
    dV_dy, dV_dx = np.gradient(potential, y, x)
    grad_x = dV_dx
    grad_y = dV_dy

    field_mag = np.hypot(ex, ey)
    grad_mag = np.hypot(grad_x, grad_y)

    denominator = np.maximum(field_mag * grad_mag, 1e-30)
    cosine = (ex * grad_x + ey * grad_y) / denominator

    valid = np.isfinite(cosine) & (field_mag > 1e-10) & (grad_mag > 1e-10)
    angle_error = np.abs(np.degrees(np.arccos(np.clip(cosine, -1, 1))) - 180.0)

    if np.any(valid):
        mean_error = float(np.mean(angle_error[valid]))
        median_error = float(np.median(angle_error[valid]))
        max_error = float(np.percentile(angle_error[valid], 95))
    else:
        mean_error = median_error = max_error = np.nan

    return angle_error, mean_error, median_error, max_error


def convergence_test(charges, xmin, xmax, ymin, ymax, resolutions):
    errors = []

    for n in resolutions:
        x = np.linspace(xmin, xmax, n)
        y = np.linspace(ymin, ymax, n)
        X, Y = np.meshgrid(x, y)

        ex, ey, potential = compute_field_potential(X, Y, charges)
        _, _, error = numerical_validation(ex, ey, potential, x, y)

        valid = np.isfinite(error)
        trimmed = error[valid]

        if trimmed.size:
            errors.append(float(np.median(trimmed)))
        else:
            errors.append(np.nan)

    return np.array(errors)


def save_field_data(x, y, ex, ey, potential):
    magnitude = np.hypot(ex, ey)

    path = OUTPUT_DIR / "field_data.csv"

    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow([
            "x_m", "y_m", "Ex_N_per_C",
            "Ey_N_per_C", "E_magnitude_N_per_C",
            "potential_V"
        ])

        for j in range(len(y)):
            for i in range(len(x)):
                writer.writerow([
                    x[i],
                    y[j],
                    ex[j, i],
                    ey[j, i],
                    magnitude[j, i],
                    potential[j, i]
                ])

    return path


def draw_charges(ax, charges):
    for xq, yq, q in charges:
        ax.scatter(
            xq,
            yq,
            s=220,
            facecolor="white",
            edgecolor="black",
            linewidth=1.5,
            zorder=10
        )

        ax.text(
            xq,
            yq,
            "+" if q > 0 else "−",
            ha="center",
            va="center",
            fontsize=14,
            fontweight="bold",
            zorder=11
        )

        ax.annotate(
            f"{q:.3e} C",
            (xq, yq),
            xytext=(8, 8),
            textcoords="offset points",
            fontsize=8
        )


def field_equipotential_plot(x, y, ex, ey, potential, charges):
    fig, ax = plt.subplots(figsize=(11, 8), dpi=160)

    finite = np.isfinite(potential)
    values = np.abs(potential[finite])
    limit = np.percentile(values, 98.5)

    masked = np.ma.masked_where(np.abs(potential) > limit, potential)

    levels = np.linspace(
        np.percentile(masked.compressed(), 2),
        np.percentile(masked.compressed(), 98),
        31
    )

    contours = ax.contour(
        x,
        y,
        masked,
        levels=levels,
        linewidths=0.8,
        alpha=0.8
    )

    ax.clabel(
        contours,
        inline=True,
        fontsize=7,
        fmt="%.2e V"
    )

    speed = np.hypot(ex, ey)

    ax.streamplot(
        x,
        y,
        ex,
        ey,
        color=np.log10(speed + 1e-30),
        density=2.0,
        linewidth=0.85,
        arrowsize=0.85
    )

    draw_charges(ax, charges)

    ax.set_title(
        "Electric Field and Equipotential Visualizer",
        fontsize=16,
        fontweight="bold",
        pad=18
    )
    ax.set_xlabel("x position (m)")
    ax.set_ylabel("y position (m)")
    ax.set_aspect("equal")
    ax.grid(alpha=0.15)

    fig.text(
        0.5,
        0.015,
        "Streamlines: electric field direction   |   Contours: constant electric potential",
        ha="center",
        fontsize=9
    )

    fig.tight_layout(rect=(0, 0.035, 1, 0.95))

    path = OUTPUT_DIR / "field_equipotential.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")

    return fig, path


def potential_plot(x, y, potential, charges):
    fig, ax = plt.subplots(figsize=(10, 8), dpi=160)

    finite = np.isfinite(potential)
    vmax = np.percentile(np.abs(potential[finite]), 98.5)
    masked = np.ma.masked_where(np.abs(potential) > vmax, potential)

    image = ax.imshow(
        masked,
        extent=[x.min(), x.max(), y.min(), y.max()],
        origin="lower",
        aspect="equal",
        cmap="coolwarm"
    )

    fig.colorbar(image, ax=ax, label="Electric potential (V)")
    draw_charges(ax, charges)

    ax.set_title(
        "Electric Potential Distribution",
        fontsize=16,
        fontweight="bold",
        pad=18
    )
    ax.set_xlabel("x position (m)")
    ax.set_ylabel("y position (m)")
    ax.set_aspect("equal")

    fig.tight_layout(rect=(0,0,1,0.95))

    path = OUTPUT_DIR / "electric_potential.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")

    return fig, path


def field_magnitude_plot(x, y, ex, ey, charges):
    magnitude = np.hypot(ex, ey)
    finite = np.isfinite(magnitude)
    vmax = np.percentile(magnitude[finite], 98.5)

    fig, ax = plt.subplots(figsize=(10, 8), dpi=160)

    masked = np.ma.masked_where(magnitude > vmax, magnitude)

    image = ax.imshow(
        masked,
        extent=[x.min(), x.max(), y.min(), y.max()],
        origin="lower",
        aspect="equal",
        norm=SymLogNorm(
            linthresh=max(vmax * 1e-5, 1e-12),
            vmax=vmax
        )
    )

    fig.colorbar(
        image,
        ax=ax,
        label="Electric field magnitude (N/C)"
    )

    draw_charges(ax, charges)

    ax.set_title(
        "Electric Field Magnitude",
        fontsize=16,
        fontweight="bold",
        pad=18
    )
    ax.set_xlabel("x position (m)")
    ax.set_ylabel("y position (m)")
    ax.set_aspect("equal")

    fig.tight_layout(rect=(0,0,1,0.95))

    path = OUTPUT_DIR / "electric_field_magnitude.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")

    return fig, path


def validation_plot(x, y, error, charges):
    positive = error[np.isfinite(error) & (error > 0)]

    if positive.size:
        lower = max(np.percentile(positive, 1), 1e-10)
        upper = max(np.percentile(positive, 99), lower * 10)
    else:
        lower, upper = 1e-10, 1.0

    masked = np.ma.masked_where(
        (~np.isfinite(error)) | (error <= 0),
        error
    )

    fig, ax = plt.subplots(figsize=(10, 8), dpi=160)

    image = ax.imshow(
        masked,
        extent=[x.min(), x.max(), y.min(), y.max()],
        origin="lower",
        aspect="equal",
        norm=SymLogNorm(
            linthresh=lower,
            vmin=lower,
            vmax=upper
        )
    )

    fig.colorbar(image, ax=ax, label="Relative error")
    draw_charges(ax, charges)

    valid = np.isfinite(error)
    median = np.median(error[valid])

    ax.text(
        0.02,
        0.98,
        f"Median relative error: {median:.3e}",
        transform=ax.transAxes,
        va="top",
        bbox=dict(
            boxstyle="round,pad=0.4",
            facecolor="white",
            alpha=0.85
        )
    )

    ax.set_title(
        r"Numerical Validation: $\mathbf{E} \approx -\nabla V$",
        fontsize=16,
        fontweight="bold",
        pad=18
    )
    ax.set_xlabel("x position (m)")
    ax.set_ylabel("y position (m)")
    ax.set_aspect("equal")

    fig.tight_layout(rect=(0,0,1,0.95))

    path = OUTPUT_DIR / "field_validation.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")

    return fig, path


def orthogonality_plot(x, y, angle_error, charges):
    finite = np.isfinite(angle_error)
    values = angle_error[finite]

    upper = max(np.percentile(values, 99), 1e-6) if values.size else 10.0

    masked = np.ma.masked_where(~finite, angle_error)

    fig, ax = plt.subplots(figsize=(10, 8), dpi=160)

    image = ax.imshow(
        masked,
        extent=[x.min(), x.max(), y.min(), y.max()],
        origin="lower",
        aspect="equal",
        vmax=upper
    )

    fig.colorbar(
        image,
        ax=ax,
        label="Deviation from 180° (degrees)"
    )

    draw_charges(ax, charges)

    ax.set_title(
        "Electric Field–Equipotential Orthogonality",
        fontsize=16,
        fontweight="bold",
        pad=18
    )
    ax.set_xlabel("x position (m)")
    ax.set_ylabel("y position (m)")
    ax.set_aspect("equal")

    fig.tight_layout(rect=(0,0,1,0.95))

    path = OUTPUT_DIR / "orthogonality_analysis.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")

    return fig, path


def convergence_plot(resolutions, errors):
    fig, ax = plt.subplots(figsize=(9, 6), dpi=160)

    valid = np.isfinite(errors) & (errors > 0)

    ax.loglog(
        np.array(resolutions)[valid],
        errors[valid],
        marker="o",
        linewidth=1.5
    )

    ax.set_title(
        "Grid-Resolution Convergence",
        fontsize=15,
        fontweight="bold"
    )
    ax.set_xlabel("Grid resolution")
    ax.set_ylabel("Median relative error")
    ax.grid(True, which="both", alpha=0.2)

    fig.tight_layout()

    path = OUTPUT_DIR / "convergence_analysis.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")

    return fig, path


def probe_analysis(charges):
    xp = read_float("Probe x position (m)")
    yp = read_float("Probe y position (m)")

    ex, ey, potential = compute_field_potential(
        np.array([[xp]]),
        np.array([[yp]]),
        charges
    )

    ex_value = float(ex[0, 0])
    ey_value = float(ey[0, 0])
    field = float(np.hypot(ex_value, ey_value))
    angle = float(np.degrees(np.arctan2(ey_value, ex_value)))
    voltage = float(potential[0, 0])

    print("\nProbe analysis")
    print(f"Position: ({xp:.6g}, {yp:.6g}) m")
    print(f"Ex: {ex_value:.6e} N/C")
    print(f"Ey: {ey_value:.6e} N/C")
    print(f"|E|: {field:.6e} N/C")
    print(f"Direction: {angle:.4f} degrees")
    print(f"Potential: {voltage:.6e} V")


def main():
    print("\n" + "=" * 66)
    print("ELECTRIC FIELD AND EQUIPOTENTIAL VISUALIZER")
    print("=" * 66)

    charges = read_charges()

    print("\nComputational domain")
    xmin = read_float("Minimum x (m)", -1.5)
    xmax = read_float("Maximum x (m)", 1.5)
    ymin = read_float("Minimum y (m)", -1.5)
    ymax = read_float("Maximum y (m)", 1.5)
    n = read_int("Grid resolution", 500, 50)

    print("\nAnalysis and output options")
    print("1  Field + equipotential map")
    print("2  Electric potential map")
    print("3  Electric field magnitude")
    print("4  Numerical validation E = -grad(V)")
    print("5  Field/equipotential orthogonality")
    print("6  Grid-resolution convergence")
    print("7  Probe-point analysis")
    print("8  Export numerical field data")
    print("9  All analyses")

    choices = {
        item.strip()
        for item in input("Enter choices separated by commas: ").split(",")
    }

    if "9" in choices:
        choices = {"1", "2", "3", "4", "5", "6", "7", "8"}

    x = np.linspace(xmin, xmax, n)
    y = np.linspace(ymin, ymax, n)
    X, Y = np.meshgrid(x, y)

    ex, ey, potential = compute_field_potential(X, Y, charges)

    figures = []

    if "1" in choices:
        figures.append(field_equipotential_plot(
            x, y, ex, ey, potential, charges
        ))

    if "2" in choices:
        figures.append(potential_plot(
            x, y, potential, charges
        ))

    if "3" in choices:
        figures.append(field_magnitude_plot(
            x, y, ex, ey, charges
        ))

    if "4" in choices:
        _, _, error = numerical_validation(
            ex, ey, potential, x, y
        )
        figures.append(validation_plot(
            x, y, error, charges
        ))

    if "5" in choices:
        angle_error, mean_error, median_error, p95_error = orthogonality_analysis(
            ex, ey, potential, x, y
        )

        print("\nOrthogonality analysis")
        print(f"Mean angular deviation: {mean_error:.6e} degrees")
        print(f"Median angular deviation: {median_error:.6e} degrees")
        print(f"95th-percentile deviation: {p95_error:.6e} degrees")

        figures.append(orthogonality_plot(
            x, y, angle_error, charges
        ))

    if "6" in choices:
        resolutions = [100, 200, 300, 400, 500, 700]
        errors = convergence_test(
            charges,
            xmin,
            xmax,
            ymin,
            ymax,
            resolutions
        )

        print("\nGrid convergence")
        for resolution, error in zip(resolutions, errors):
            print(f"{resolution:4d} × {resolution:<4d} : {error:.6e}")

        figures.append(convergence_plot(
            resolutions,
            errors
        ))

    if "7" in choices:
        probe_analysis(charges)

    if "8" in choices:
        path = save_field_data(
            x, y, ex, ey, potential
        )
        print(f"\nSaved numerical data: {path}")

    if figures:
        print("\nGenerated figures:")
        for _, path in figures:
            print(path)

        plt.show()


if __name__ == "__main__":
    main()
