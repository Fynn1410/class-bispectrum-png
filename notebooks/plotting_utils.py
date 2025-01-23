import numpy as np
import matplotlib.pyplot as plt

# This function is used to "cut" a function in regions of positive and negative values
# such that we can plot it in log log
def plot_segments(ax, x, y, color, legend, ls, use_labels):
    label_counter=0
    segments = []
    current_segment = [0]  # Start from the first index

    for i in range(1, len(y)):
        if (y[i] > 0) != (y[i - 1] > 0):  # Sign change
            current_segment.append(i)
            segments.append(current_segment)
            current_segment = [i]
    current_segment.append(len(y))  # Add the last segment
    segments.append(current_segment)

    for segment in segments:
        start, end = segment[0], segment[1]
        x_segment = x[start:end]
        y_segment = y[start:end]
        if np.all(y_segment > 0):
            if use_labels and label_counter<1:
                ax.plot(x_segment, y_segment, linestyle=ls, color=color, label=rf"{legend}")
            else:
                ax.plot(x_segment, y_segment, linestyle=ls, color=color)
            label_counter+=1
        elif np.all(y_segment < 0):
            ax.plot(x_segment, np.abs(y_segment), linestyle='--', color=color)
    
    if label_counter==0:
        print(rf"You don't use any labels for {legend}")



def plot_compare(x, y11, y12, y21, y22, y31, y32, xlabel,ylabel, label1, label2, titles, save_name, colors, lin_log, abs_rel, sharexy):
    # Create a 2x3 grid of subplots with shared axes
    if sharexy:
        fig, ax = plt.subplots(2, 3, figsize=(12, 16/3), sharex=True, sharey='row', gridspec_kw={'height_ratios': [3, 1], 'hspace': 0, 'wspace': 0})
    else: 
        fig, ax = plt.subplots(2, 3, figsize=(12, 16/3), gridspec_kw={'height_ratios': [3, 1], 'hspace': 0, 'wspace': 0})

    res_y_lim = max([max(np.abs(y11-y12)/y12), max(np.abs(y21-y22)/y22), max(np.abs(y31-y32)/y32)])*200
    # Top plots: the two curves
    
    for i, (y1, y2) in enumerate([(y11, y12), (y21, y22), (y31, y32)]):
        if i==0:
            ax[0, i].plot(x, y1, linewidth=2, color=colors[0], ls="solid", label=label1)
            ax[0, i].plot(x, y2, linewidth=1, color=colors[1], ls="dashed", label=label2)
        else:
            ax[0, i].plot(x, y1, linewidth=2, color=colors[0], ls="solid")
            ax[0, i].plot(x, y2, linewidth=1, color=colors[1], ls="dashed")
        #ax[0, i].set_yscale("log")
        ax[0, i].set_title(f"{titles[i]}")  # Add unique titles for subplots
        if lin_log[0] == "log":
                ax[1, i].set_xscale("log")
        if lin_log[1] == "log":
                ax[0, i].set_yscale("log")

    # Shared y-label for the upper row
    ax[0, 0].set_ylabel(ylabel)
    #ax[0, 0].legend(loc="best")
    

    # Bottom plots: residuals
    for i, (y1, y2) in enumerate([(y11, y12), (y21, y22), (y31, y32)]):
        if abs_rel=="abs" or abs_rel=="absolute":
            residual = y2 - y1
        elif abs_rel =="rel" or abs_rel=="relative":
            residual = 100 * (y2 - y1) / y2
        else:
            residual=0
            print("abs_rel is invalid:", residual)
        ax[1, i].plot(x, residual, color="black")
        ax[1, i].axhline(0, color="black", linestyle="--", linewidth=1)
        ax[1, i].set_xlabel(xlabel)
        ax[1, i].grid(axis="y")
    ax[1, 0].set_ylim(-res_y_lim, res_y_lim)


    # Shared y-label for the bottom row
    if abs_rel =="rel" or abs_rel=="relative": 
        ax[1, 0].set_ylabel(rf"Res. $\left[  \%\right]$")
    elif abs_rel=="abs" or abs_rel=="absolute":
        ax[1, 0].set_ylabel(rf"Abs. Res.")
    
    #ax[1, 0].set_ylim(-1.5,1.5)
    fig.legend(loc='center', bbox_to_anchor=(0.5, -0.015), ncol=3)

    # Tight layout for better spacing
    # plt.tight_layout()

    # Save and show the figure
    plt.savefig(f"Dennis_Fynn_Comparison/{save_name}.pdf", bbox_inches='tight')
    plt.show()
