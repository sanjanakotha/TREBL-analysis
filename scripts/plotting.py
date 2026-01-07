import matplotlib.pyplot as plt
import seaborn as sns

def plot_reads_histogram(map_df, save_path=None, ax=None, **kwargs):    
    """
    Plots a histogram of read counts on a log-log scale.

    Parameters
    ----------
    map_df : pd.DataFrame
        DataFrame containing at least a "count" column.
    save_path : str, optional
        Path to save the figure.
    ax : matplotlib.axes.Axes, optional
        Axis to plot on. If None, a new figure is created.
    **kwargs : dict
        Additional keyword arguments passed to sns.histplot.
    
    Returns
    -------
    matplotlib.axes.Axes
        Axis with the histogram plotted.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4), dpi=300)
    
    sns.histplot(map_df["count"], log_scale=(True, True), ax=ax, **kwargs)
    
    # Ensure log scale is applied
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel("Read Count")
    
    sns.despine()
    
    if save_path:
        ax.get_figure().savefig(save_path, bbox_inches="tight")
    
    return ax
    
    
    

