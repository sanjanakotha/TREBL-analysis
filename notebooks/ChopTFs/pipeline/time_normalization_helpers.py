import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import tqdm
from scipy.stats import norm, kstest

def plot_hists_over_time(activities, column = "activity_directional", sample=False, kde=False, xlabel = "Activity", stat = 'density'):
    sns.set_context('talk')

    times = np.array(activities["time"].drop_duplicates().sort_values())
    
    fig,axs = plt.subplots(len(times),1, sharex = True, figsize = (12,12), dpi = 300)
    
    for i in tqdm.tqdm(range(len(times))):
        df = activities[activities["time"] == times[i]]

        if sample:
            df = df.sample(sample)

        if not kde:
            sns.histplot(data = df, x = column, ax = axs[i], stat = stat, label = "All Tiles")
        else:
            sns.kdeplot(data = df, x = column, ax = axs[i])
            
        axs[i].set_xlabel("")
        axs[i].set_ylabel(f"{times[i]} min.", rotation = 0, va = 'center', ha = 'right', labelpad = 15)
        
    fig.supxlabel(xlabel, y = 0.05)
    fig.supylabel("Count", x = -0.07)
    sns.despine()
    return fig,axs

def fit_inactive_gaussian_progressive_fast(df, ax=None):
    scores = np.array(df["score"])
    search_ceiling = np.percentile(scores, 99)
    truncated = scores[scores < search_ceiling]

    s = np.sort(truncated)
    n_total = len(s)

    # cumulative stats
    cumsum = np.cumsum(s)
    cumsum_sq = np.cumsum(s**2)
    n = np.arange(1, n_total + 1)

    means = cumsum / n
    vars_ = (cumsum_sq / n) - means**2
    stds = np.sqrt(np.maximum(vars_, 1e-12))

    best_score = np.inf
    best_mu = None

    min_points = max(50, int(n_total * 0.02))
    step = max(10, n_total // 200)

    worse_count = 0

    for i in range(min_points, n_total, step):
        mu = means[i-1]
        sigma = stds[i-1]

        if sigma == 0:
            continue

        window = s[:i]

        # FAST proxy instead of KS
        z = (window - mu) / sigma
        skew = np.mean(z**3)
        kurt = np.mean(z**4)

        score = abs(skew) + abs(kurt - 3)

        if score < best_score:
            best_score = score
            best_mu = mu
            worse_count = 0
        else:
            worse_count += 1

        if worse_count > 20:
            break

    peak_location = best_mu

    left_of_peak = s[s < peak_location]
    mirrored = np.concatenate([left_of_peak, 2 * peak_location - left_of_peak])

    mu = peak_location
    sigma = np.std(mirrored)

    if ax is not None:
        sns.histplot(truncated, ax=ax, stat="density", bins=80, kde=True)
        ax.axvline(mu, color="red", linestyle="dashed", lw=2)

        x = np.linspace(min(truncated), max(truncated), 500)
        ax.plot(x, norm.pdf(x, mu, sigma), "r-", lw=2)

    df = df.copy()
    df["z-score"] = (df["score"] - mu) / sigma

    return mu, sigma, truncated, df

def plot_ridge(activities, column="activity_directional", sample=False, bw_adjust=1.0, xlabel = "Activity"):
    sns.set_context('talk')
    times = np.array(activities["time"].drop_duplicates().sort_values())  # reverse: earliest on bottom
    n = len(times)

    palette = sns.color_palette("crest", n)
    fig, axes = plt.subplots(n, 1, figsize=(10, n * 1.2), sharex = True, sharey = True, dpi = 300)
    axes = np.atleast_1d(axes)

    x_min = activities[column].quantile(0.01)
    x_max = activities[column].quantile(0.99)
    bins = np.linspace(x_min, x_max, 100)

    for i, (ax, t) in enumerate(zip(axes, times)):
        df = activities[activities["time"] == t]
        if sample:
            df = df.sample(sample)

        data = df[column].dropna()
        ax.hist(data, bins=bins, color=palette[i], alpha=1, histtype='stepfilled', linewidth=1.5, edgecolor='white', density=True)
        ax.set_xlim(x_min, x_max)
        ax.set_ylabel(f"{t} min.", rotation=0, va='bottom', ha='right', labelpad=15, fontsize='medium', y = 0)
        ax.set_xlabel("")
        ax.set_yticks([])
        ax.set_facecolor("none")
        sns.despine(ax=ax, left=True)

    fig.supxlabel(xlabel, y=0.02)
    plt.subplots_adjust(hspace=-0.4)

    return fig, axes