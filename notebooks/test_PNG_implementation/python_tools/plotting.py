import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec








def plot_comparison(n_triag,
                    y11, y12, y21, y22, y31, y32,
                    k1_arr, k2_arr, k3_arr,
                    h,
                    label1, label2):
    """
    y11 = B0_model
    y21 = B2_model
    y31 = B4_model

    y12 = B0_reference
    y22 = B2_reference
    y32 = B4_reference

    Input units:
    k: 1/Mpc
    B: Mpc^6

    Plot Units:
    k: h/Mpc
    B: (Mpc/h)^6
    """

    colors = [
        "#332288",  # deep purple
        "#29A3E1",  # vibrant blue (your choice)
        "#238B45",  # fresher green (darker, cooler than #117733)
        "#FDAE61",  # bright orange (replaces #FDB863 for more vibrance)
        "#D73027",  # vivid red (stronger separation from green)
        "#1F9FB6",  # refined aqua blue (replaces cyan, brighter and cleaner)
    ]



    colors_B = [
        colors[0],
        colors[1]
    ]


    colors_k = [
        "tab:blue",
        "tab:red",
        "tab:green"
    ]

    titles = ["Monopole", "Quadrupole", "Hexadecapole"]
    ylabel = rf"$B_\ell(k_1, k_2, k_3)$ $\left[(h/\mathrm{{Mpc}})^{{-6}}\right]$"
    xlabel = "Triangle Counter"


    # Settings
    res_y_lim = max([
        max(np.abs(y11 - y12) / y12),
        max(np.abs(y21 - y22) / y22),
        max(np.abs(y31 - y32) / y32)
    ]) * 200

    # === Create figure with 4-row GridSpec (Top, Residuals, Spacer, k-plot) ===
    fig = plt.figure(figsize=(9, 5.25))
    gs = GridSpec(4, 3, figure=fig, height_ratios=[3, 1, 0.4, 1], hspace=0.0, wspace=0.0)

    # === Add top and bottom row axes (Plots 1–6) ===
    axes_top = []
    axes_bottom = []

    for col in range(3):
        if col == 0:
            ax_top = fig.add_subplot(gs[0, col])
            ax_bottom = fig.add_subplot(gs[1, col], sharex=ax_top)
        else:
            ax_top = fig.add_subplot(gs[0, col], sharex=axes_top[0], sharey=axes_top[0])
            ax_bottom = fig.add_subplot(gs[1, col], sharex=axes_top[0])
        axes_top.append(ax_top)
        axes_bottom.append(ax_bottom)

    # === Add Plot 7 (full-width) ===
    ax_k = fig.add_subplot(gs[3, :]) # sharex=axes_top[0]

    # === Plot Plots 1–6 ===
    for i, (ax_top, ax_bot, y1, y2, title) in enumerate(zip(
        axes_top, axes_bottom,
        [y11, y21, y31], [y12, y22, y32], titles
    )):
        # Top: main comparison
        ax_top.plot(n_triag, y1 * h**6, color=colors_B[0], label=label1 if i == 0 else "")
        ax_top.plot(n_triag, y2 * h**6, color=colors_B[1], ls="--", label=label2 if i == 0 else "")
        ax_top.set_yscale("log")
        ax_top.set_title(title)

        if i == 0:
            ax_top.set_ylabel(ylabel)
        else:
            ax_top.tick_params(labelleft=False)

        # Bottom: residuals
        residual = 100 * (y2 - y1) / y1
        print(f"residual max deviation: {np.max(np.abs(residual)):.3f}%")
        ax_bot.axhline(0, color="black", linestyle="--", linewidth=1)
        ax_bot.plot(n_triag, residual, color="black")
        ax_bot.set_ylim(-res_y_lim, res_y_lim)
        ax_bot.grid(axis="y")
        ax_bot.set_xlabel(xlabel)
        ax_bot.set_xticks([0, 2000, 4000])
        # ax_bot.set_yticks([-0.1, 0, 0.1])

        if i == 0:
            ax_bot.set_ylabel(rf"Res. $\left[  \%\right]$")
        else:
            ax_bot.tick_params(labelleft=False)

    # === Plot 7 ===
    ax_k.plot(n_triag, k3_arr / h, color=colors_k[2], label="k3")
    ax_k.plot(n_triag, k2_arr / h, color=colors_k[1], label="k2")
    ax_k.plot(n_triag, k1_arr / h, color=colors_k[0], label="k1")
    ax_k.set_ylabel(r"$k \, [h/\mathrm{Mpc}]$")
    ax_k.set_xlabel("Triangle Counter")

    ax_k.set_xticks([0, 1000, 2000, 3000, 4000, 5000])


    # === Final touches ===
    fig.legend(loc='center', bbox_to_anchor=(0.54, 0.0031), ncol=3)
    fig.tight_layout(rect=[0, 0.03, 1, 1])  # Keep space for legend
    plt.show()


def plot_multipoles(n_triag,
                    y1, y2, y3,
                    k1_arr, k2_arr, k3_arr,
                    h,
                    label,
                    global_title=None):
    """
    y1 = B0
    y2 = B2
    y3 = B4

    y1, y2, y3 may each be a single array OR a list/tuple of arrays (one per
    curve to overlay). `label` may be a single string or a list of strings,
    one per curve.

    Input units:
    k: 1/Mpc
    B: Mpc^6

    Plot Units:
    k: h/Mpc
    B: (Mpc/h)^6
    """

    colors = [
        "#332288",  # deep purple
        "#29A3E1",  # vibrant blue
        "#238B45",  # fresher green
        "#FDAE61",  # bright orange
        "#D73027",  # vivid red
        "#1F9FB6",  # refined aqua blue
    ]

    colors_k = [
        "tab:blue",
        "tab:red",
        "tab:green"
    ]

    titles = ["Monopole", "Quadrupole", "Hexadecapole"]
    ylabel = rf"$B_\ell(k_1, k_2, k_3)$ $\left[(h/\mathrm{{Mpc}})^{{-6}}\right]$"
    xlabel = "Triangle Counter"

    def _to_list(y):
        return list(y) if isinstance(y, (list, tuple)) else [y]

    y1_list = _to_list(y1)
    y2_list = _to_list(y2)
    y3_list = _to_list(y3)
    n_curves = len(y1_list)
    assert len(y2_list) == n_curves and len(y3_list) == n_curves, \
        "y1, y2, y3 must contain the same number of arrays"

    labels = list(label) if isinstance(label, (list, tuple)) else [label]
    labels += [""] * (n_curves - len(labels))

    # === Create figure with 3-row GridSpec (Top, Spacer, k-plot) ===
    fig = plt.figure(figsize=(9, 5.5))
    gs = GridSpec(3, 3, figure=fig, height_ratios=[3, 0.7, 1], hspace=0.0, wspace=0.0)

    # === Add top row axes (Plots 1–3) ===
    axes_top = []

    for col in range(3):
        if col == 0:
            ax_top = fig.add_subplot(gs[0, col])
        else:
            ax_top = fig.add_subplot(gs[0, col], sharex=axes_top[0], sharey=axes_top[0])
        axes_top.append(ax_top)

    # === Add k-plot (full-width) ===
    ax_k = fig.add_subplot(gs[2, :])

    # === Plot multipoles ===
    any_negative = any(
        np.any(np.asarray(y) * h**6 < 0)
        for ys in [y1_list, y2_list, y3_list] for y in ys
    )
    n_arr = np.asarray(n_triag)
    stride = max(1, len(n_arr) // 50)

    for i, (ax_top, ys, title) in enumerate(zip(
        axes_top, [y1_list, y2_list, y3_list], titles
    )):
        for j, y in enumerate(ys):
            y_scaled = np.asarray(y) * h**6
            color = colors[j % len(colors)]
            curve_label = labels[j] if i == 0 else ""
            ax_top.plot(n_arr, np.abs(y_scaled), color=color, label=curve_label)

            neg_idx = np.where(y_scaled < 0)[0][::stride]
            ax_top.plot(n_arr[neg_idx], np.abs(y_scaled)[neg_idx],
                        linestyle="none", marker="x", markersize=4, color=color)

        ax_top.set_yscale("log")
        ax_top.set_title(title)
        ax_top.set_xticks([0, 2000, 4000])

        if i == 0:
            ax_top.set_ylabel(ylabel)
            if any_negative:
                ax_top.plot([], [], linestyle="none", marker="x", markersize=4,
                            color="black", label="negative")
        else:
            ax_top.tick_params(labelleft=False)

    # === k-plot ===
    ax_k.plot(n_triag, k3_arr / h, color=colors_k[2], label="k3")
    ax_k.plot(n_triag, k2_arr / h, color=colors_k[1], label="k2")
    ax_k.plot(n_triag, k1_arr / h, color=colors_k[0], label="k1")
    ax_k.set_yscale("log")
    ax_k.set_ylabel(r"$k \, [h/\mathrm{Mpc}]$")
    ax_k.set_xlabel(xlabel)

    ax_k.set_xticks([0, 1000, 2000, 3000, 4000, 5000])

    # === Final touches ===
    if global_title is not None:
        fig.suptitle(global_title, fontsize='medium', y=0.98)
    fig.legend(loc='center', bbox_to_anchor=(0.54, 0.0031), ncol=4)
    top_rect = 0.97 if global_title is not None else 1.0
    fig.tight_layout(rect=[0, 0.07, 1, top_rect])
    plt.show()

