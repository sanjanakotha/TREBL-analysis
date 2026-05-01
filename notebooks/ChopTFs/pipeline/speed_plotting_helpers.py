import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns

sns.set_style('ticks')

def plot_model_fits(
    activities,
    fit_results,
    sample_adseqs,
    model_fn,
    model_params,
    y_col="activity",
    yerr_col=None,
    ci_low_col="ci_low",
    ci_hi_col="ci_hi",
    title="Model Fits",
    n_cols=2,
    n_rows=5,
    figsize =(6,8),
    data_color="steelblue",
    fit_color="tomato",
):

    if len(sample_adseqs) > 10:
        sns.set_context("talk")
    
    
        n = len(sample_adseqs)
        t_fine = np.linspace(activities["time"].min(), activities["time"].max(), 200)
    
        fig, axes = plt.subplots(
            n_rows, n_cols,
            figsize=figsize,
            sharex=True,
            sharey=True,
            dpi=300
        )
        axes = np.array(axes).flatten()
    
        for ax, adseq in zip(axes, sample_adseqs):
            grp = activities[activities["ADseq"] == adseq].sort_values("time")
            t = grp["time"].values
            y = grp[y_col].values
    
            # --- error bars ---
            if yerr_col is not None:
                yerr = grp[yerr_col].values
            else:
                if ci_low_col is None or ci_hi_col is None:
                    raise ValueError("Could not find CI columns. Pass ci_low_col and ci_hi_col.")
    
                yerr = [
                    y - grp[ci_low_col].values,
                    grp[ci_hi_col].values - y
                ]
    
            ax.errorbar(
                t, y,
                yerr=yerr,
                fmt="o",
                color=data_color,
                ecolor=data_color,
                elinewidth=1.5,
                capsize=2,
                markersize=4,
                alpha=0.9
            )
    
            # --- model fit ---
            row = fit_results[fit_results["ADseq"] == adseq]
            if len(row) > 0:
                params = [row.iloc[0][p] for p in model_params]
                ax.plot(
                    t_fine,
                    model_fn(t_fine, *params),
                    color=fit_color,
                    linewidth=2,
                )
    
        # hide unused axes
        for ax in axes[n:]:
            ax.set_visible(False)
    
        sns.despine()
        plt.tight_layout()
    
        data_handle = mlines.Line2D([], [], color=data_color, marker="o",
                                    linestyle="None", label="data ± CI")
        fit_handle = mlines.Line2D([], [], color=fit_color, linewidth=2, label="fit")
    
        # fig.legend(handles=[data_handle, fit_handle],
        #            loc="upper right", bbox_to_anchor=(1, 1.07), ncols = 2)
    
        plt.suptitle(title, y=1.05)
        fig.supxlabel("Time", y=-0.05)
        fig.supylabel("Activity", x=-0.03)
    
        plt.show()