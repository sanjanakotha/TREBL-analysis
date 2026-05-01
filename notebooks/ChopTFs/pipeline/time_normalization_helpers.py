import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import tqdm
from scipy.stats import norm, kstest
import pandas as pd

def plot_hists_over_time(activities, column = "activity_directional", sample=False, kde=False, xlabel = "Activity", stat = 'density',step = True, xlims = (0,100)):
    sns.set_context('talk')

    times = np.array(activities["time"].drop_duplicates().sort_values())
    fig,axs = plt.subplots(len(times),1, sharex = True, sharey = True,figsize = (8,12), dpi = 300)

    lower, upper = np.percentile(activities[column].dropna(), xlims)
    
    for i in tqdm.tqdm(range(len(times))):
        df = activities[activities["time"] == times[i]]

        if sample:
            df = df.sample(sample)

        if not kde:
            if step:
                sns.histplot(data = df, x = column, ax = axs[i], stat = stat, label = "All Tiles", element = "step", alpha = 0.5)
            else:
                sns.histplot(data = df, x = column, ax = axs[i], stat = stat, label = "All Tiles")
        else:
            sns.kdeplot(data = df, x = column, ax = axs[i])
            
        axs[i].set_xlabel("")
        axs[i].set_ylabel(f"{times[i]} min.", rotation = 0, va = 'center', ha = 'right', labelpad = 15)

        #axs[i].set_xlim(lower, upper)
        
    fig.supxlabel(xlabel, y = 0.05)
    fig.supylabel("Count", x = -0.1)
    sns.despine()
    return fig,axs

def plot_gaussian_fits_over_time(activities, cached_fits, column = "activity_directional", xlabel = "Activity", dist_color = 'C0', xlims = (0,100)):
    sns.set_context('talk')

    times = np.array(activities["time"].drop_duplicates().sort_values())
    fig,axs = plt.subplots(len(times),1, sharex = True, figsize = (8,12), dpi = 300, sharey = True)

    lower, upper = np.percentile(activities[column].dropna(), xlims)
    
    for i in tqdm.tqdm(range(len(times))):
        time = times[i]
        df = activities[activities["time"] == time]

        # get cached fit for this (rep, time)
        current_rep = df["rep"].iloc[0]
        mu, sigma, truncated, _ = cached_fits[(current_rep, time)]

        # histogram (same style as before)
        sns.histplot(
            data=df,
            x=column,
            ax=axs[i],
            stat='density',
            label="All Tiles",
            element="step",
            alpha=0.5,
            color = dist_color
        )

        # gaussian fit
        x = np.linspace(lower, upper, 200)

        from scipy.stats import norm
        
        cdf_upper = norm.cdf(upper, mu, sigma)
        cdf_lower = norm.cdf(lower, mu, sigma)
        Z = cdf_upper - cdf_lower  # normalization constant
        
        axs[i].plot(x, norm.pdf(x, mu, sigma) / Z, color='black', lw=2)

        #axs[i].plot(x, norm.pdf(x, mu, sigma), color='red', lw=2)
        axs[i].axvline(mu, color='black', linestyle='dashed', lw=2)

        axs[i].set_xlabel("")
        axs[i].set_ylabel(
            f"{time} min.",
            rotation=0,
            va='center',
            ha='right',
            labelpad=15
        )

        #axs[i].set_xlim(lower, upper)

    
    fig.supxlabel(xlabel, y = 0.05)
    fig.supylabel("Count", x = -0.1)
    sns.despine()
        
    return fig,axs

def fit_inactive_gaussian_progressive_fast(df, plot_steps=False):
    scores = np.array(df["score"])
    mu, sigma = None, None
    
    # --- STEP 1: raw scores ---
    if plot_steps:
        plt.figure(figsize=(10,5), dpi=150)
        sns.histplot(scores, bins=80, stat='density', color='skyblue')
        plt.title("1. Raw scores")
        plt.xlabel("Score")
        xlims = plt.xlim()  # save x-limits to reuse
        plt.show()
    else:
        xlims = (np.min(scores), np.max(scores))
    
    # --- STEP 2: truncate top 1% ---
    search_ceiling = np.percentile(scores, 99)
    truncated = scores[scores < search_ceiling]
    
    if plot_steps:
        plt.figure(figsize=(10,5), dpi=150)
        sns.histplot(truncated, bins=80, stat='density', color='lightgreen')
        plt.axvline(search_ceiling, color='red', linestyle='--', label='99th percentile')
        plt.title("2. Truncated scores (remove top 1%)")
        plt.xlim(xlims)
        plt.show()
    
    # --- STEP 3: progressive Gaussian fits ---
    s = np.sort(truncated)
    n_total = len(s)
    
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

    scores_trace = []
    mus_trace = []
    
    if plot_steps:
        plt.figure(figsize=(10,5), dpi=150)
        counts, bins, _ = plt.hist(s, bins=80, density=True, color='lightgray', alpha=0.5, label='Truncated data')
    
    for i in range(min_points, n_total, step):
        mu_i = means[i-1]
        sigma_i = stds[i-1]
        if sigma_i == 0:
            continue
    
        window = s[:i]
        z = (window - mu_i) / sigma_i
        skew = np.mean(z**3)
        kurt = np.mean(z**4)
        score = abs(skew) + abs(kurt - 3)
    
        scores_trace.append(score)
        mus_trace.append(mu_i)
    
        # overlay example fits
        if plot_steps and i % (5*step) == 0:
            x_fit = np.linspace(bins[0], bins[-1], 200)       # use histogram bins
            y_fit = norm.pdf(x_fit, mu_i, sigma_i)
            y_fit *= max(counts) / max(y_fit)                 # scale to histogram
            plt.plot(x_fit, y_fit, color='red', alpha=0.4)
    
        if score < best_score:
            best_score = score
            best_mu = mu_i
            worse_count = 0
        else:
            worse_count += 1
        if worse_count > 100000:
            break
    
    if plot_steps:
        plt.title("3. Progressive Gaussian fits (example overlaid)")
        plt.xlabel("Score")
        plt.xlim(xlims)
        plt.show()
    peak_location = best_mu

    # --- STEP 4: mirror left of peak ---
    left_of_peak = s[s < peak_location]
    mirrored = np.concatenate([left_of_peak, 2 * peak_location - left_of_peak])
    mu = peak_location
    sigma = np.std(mirrored)
    
    if plot_steps:
        plt.figure(figsize=(10,5), dpi=150)
        sns.histplot(mirrored, bins=80, stat='density', color='orange')
        x = np.linspace(xlims[0], xlims[1], 200)
        plt.plot(x, norm.pdf(x, mu, sigma), 'r-', label='Mirrored Gaussian')
        plt.title("4. Fit normal distribution with data mirrored to left of peak")
        plt.xlabel("Score")
        plt.xlim(xlims)
        plt.show()
    
    # --- STEP 5: final z-score ---
    df = df.copy()
    df["z-scored_activity"] = (df["score"] - mu) / sigma
    df["shifted_activity"] = df["score"] - mu
    df["mu"] = mu
    df["sigma"] = sigma
    
    if plot_steps:
        plt.figure(figsize=(10,5), dpi=150)
        sns.histplot(df["z-scored_activity"], bins=80, stat='density', color='purple')
        plt.title("5. Final z-scored activity")
        plt.xlabel("Z-score")
        plt.xlim((xlims - mu) / sigma)
        plt.show()
    
    return mu, sigma, truncated, df

def fit_inactive_gaussian_progressive_fast_right_anchored(df, plot_steps=False):
    scores = np.array(df["score"])

    # --- STEP 1: raw scores ---
    if plot_steps:
        plt.figure(figsize=(10,5), dpi=150)
        sns.histplot(scores, bins=80, stat='density', color='skyblue')
        plt.title("1. Raw scores")
        xlims = plt.xlim()
        plt.show()
    else:
        xlims = (np.min(scores), np.max(scores))

    # --- STEP 2: truncate top 1% ---
    search_ceiling = np.percentile(scores, 99)
    truncated = scores[scores < search_ceiling]

    if plot_steps:
        plt.figure(figsize=(10,5), dpi=150)
        sns.histplot(truncated, bins=80, stat='density', color='lightgreen')
        plt.axvline(search_ceiling, color='red', linestyle='--', label='99th percentile')
        plt.title("2. Truncated scores (remove top 1%)")
        plt.xlim(xlims)
        plt.show()

    # --- STEP 3: progressive fits anchored at RIGHT end, expanding left ---
    s = np.sort(truncated)[::-1]   # sort descending so s[0] is the top, s[:i] expands leftward
    n_total = len(s)

    best_score = np.inf
    best_idx = None
    min_points = max(50, int(n_total * 0.02))
    step = max(10, n_total // 200)
    worse_count = 0

    # precompute running stats from the right
    cumsum = np.cumsum(s)
    cumsum_sq = np.cumsum(s**2)
    n = np.arange(1, n_total + 1)
    means = cumsum / n
    vars_ = (cumsum_sq / n) - means**2
    stds = np.sqrt(np.maximum(vars_, 1e-12))

    if plot_steps:
        plt.figure(figsize=(10,5), dpi=150)
        counts, bins, _ = plt.hist(truncated, bins=80, density=True, color='lightgray', alpha=0.5)

    for i in range(min_points, n_total, step):
        mu_i = means[i-1]
        sigma_i = stds[i-1]
        if sigma_i == 0:
            continue

        window = s[:i]
        z = (window - mu_i) / sigma_i
        skew = np.mean(z**3)
        kurt = np.mean(z**4)
        score = abs(skew) + abs(kurt - 3)

        if plot_steps and i % (5*step) == 0:
            x_fit = np.linspace(bins[0], bins[-1], 200)
            y_fit = norm.pdf(x_fit, mu_i, sigma_i)
            y_fit *= max(counts) / max(y_fit)
            plt.plot(x_fit, y_fit, color='red', alpha=0.4)

        if score < best_score:
            best_score = score
            best_idx = i
            worse_count = 0
        else:
            worse_count += 1
        if worse_count > 100000:
            break

    if plot_steps:
        plt.title("3. Progressive Gaussian fits anchored at right, expanding left")
        plt.xlabel("Score")
        plt.xlim(xlims)
        plt.show()

    # best window is s[:best_idx] (right-anchored, expanding left)
    best_window = s[:best_idx]
    peak_location = means[best_idx - 1]

    # --- STEP 4: mirror right of peak (since we're anchored on the right) ---
    right_of_peak = best_window[best_window > peak_location]
    mirrored = np.concatenate([right_of_peak, 2 * peak_location - right_of_peak])
    mu = peak_location
    sigma = np.std(mirrored)

    if plot_steps:
        plt.figure(figsize=(10,5), dpi=150)
        sns.histplot(mirrored, bins=80, stat='density', color='orange')
        x = np.linspace(xlims[0], xlims[1], 200)
        plt.plot(x, norm.pdf(x, mu, sigma), 'r-', label='Mirrored Gaussian')
        plt.title("4. Fit normal distribution with data mirrored to right of peak")
        plt.xlabel("Score")
        plt.xlim(xlims)
        plt.show()

    # --- STEP 5: final z-score ---
    df = df.copy()
    df["z-scored_activity"] = (df["score"] - mu) / sigma
    df["shifted_activity"] = df["score"] - mu
    df["mu"] = mu
    df["sigma"] = sigma

    if plot_steps:
        plt.figure(figsize=(10,5), dpi=150)
        sns.histplot(df["z-scored_activity"], bins=80, stat='density', color='purple')
        plt.title("5. Final z-scored activity")
        plt.xlabel("Z-score")
        plt.show()

    return mu, sigma, truncated, df


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import norm

def fit_normal_dbd(df, Q=0.3, plot_steps=False):
    scores = np.array(df["score"])
    xlims = (np.min(scores), np.max(scores))

    # --- STEP 1: raw scores ---
    if plot_steps:
        plt.figure(figsize=(10, 5), dpi=150)
        sns.histplot(scores, bins=80, stat='density', color='skyblue')
        plt.title("1. Raw scores")
        plt.show()

    # --- STEP 2: iterative trimming (ROUT) with KDE initialization ---
    working = scores.copy()
    for i in range(100):
        if i == 0:
            kde = stats.gaussian_kde(working, bw_method='silverman')
            x_grid = np.linspace(working.min(), working.max(), 1000)
            mu_iter = x_grid[np.argmax(kde(x_grid))]
            # use IQR-based sigma instead of std to avoid tail inflation
            iqr = np.percentile(working, 75) - np.percentile(working, 25)
            sigma_iter = iqr / 1.35

        else:
            mu_iter, sigma_iter = stats.norm.fit(working)

        n = len(working)
        p_values = 2 * stats.norm.sf(np.abs(working - mu_iter) / sigma_iter)
        sorted_idx = np.argsort(p_values)
        sorted_p = p_values[sorted_idx]
        bh_threshold = (np.arange(1, n + 1) / n) * Q
        passes = sorted_p <= bh_threshold
        if not passes.any():
            break
        last_outlier = np.where(passes)[0][-1]
        keep_mask = np.ones(n, dtype=bool)
        keep_mask[sorted_idx[:last_outlier + 1]] = False
        if keep_mask.all():
            break
        working = working[keep_mask]

    # after ROUT converges and you have mu
    mu, sigma = stats.norm.fit(working)

    if plot_steps:
        plt.figure(figsize=(10, 5), dpi=150)
        sns.histplot(scores, bins=80, stat='density', color='lightgray', alpha=0.6, label='all data')
        sns.histplot(working, bins=80, stat='density', color='salmon', alpha=0.6, label='inliers')
        x = np.linspace(xlims[0], xlims[1], 400)
        plt.plot(x, norm.pdf(x, mu, sigma), 'r-', lw=2, label=f'fit: μ={mu:.3f}, σ={sigma:.3f}')
        plt.axvline(mu, color='darkred', linestyle='--')
        plt.title(f"2. ROUT iterative trimming (Q={Q})")
        plt.xlim(xlims)
        plt.legend()
        plt.show()

    # --- STEP 3: final z-score ---
    df = df.copy()
    df["z-scored_activity"] = (df["score"] - mu) / sigma
    df["shifted_activity"]  = df["score"] - mu
    df["mu"]                = mu
    df["sigma"]             = sigma
    if plot_steps:
        plt.figure(figsize=(10, 5), dpi=150)
        sns.histplot(df["z-scored_activity"], bins=80, stat='density', color='purple')
        plt.title("3. Final z-scored activity")
        plt.xlabel("Z-score")
        plt.show()

    return mu, sigma, working, df
    
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

def plot_all_control_groups(all_activities, reps = [2,3,4], max_activity = np.inf):
    chopTF_controls = pd.read_csv("../../../data/chopTFs_controls.csv")
    all_activities = all_activities.rename(columns = {"seq" : "ADseq"})
    all_activities = all_activities[all_activities["time"] <= max_activity]


    control_activities = pd.merge(all_activities, chopTF_controls[["ADseq", "Name", "group"]], on = "ADseq")

    control_activities_grouped = control_activities.groupby(["ADseq", "rep", "time", "Name", "group"])[["z-scored_activity", "shifted_activity"]].agg(["mean", "std"]).reset_index()
    control_activities_grouped.columns = ['_'.join(col).strip() if isinstance(col, tuple) else col for col in control_activities_grouped.columns]


    for group_name in control_activities_grouped.sort_values(by = "group_")["group_"].unique():
    
        
        one_group_activities = control_activities_grouped[control_activities_grouped["group_"] == group_name]
        
        # Ensure your dataframe is sorted by time
        df_plot = one_group_activities.sort_values("time_")  # replace 'df' with your dataframe
        
        n_rows, n_cols = 2, 3
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(12,6), sharex=True, sharey='row', dpi = 300)
        
        # Generate a color palette for ADseqs
        adseqs = df_plot["ADseq_"].unique()
        palette = sns.color_palette("tab10", n_colors=len(adseqs))
        adseq_color = dict(zip(adseqs, palette))
        
        # Map ADseq_ → Name_ for legend labels
        adseq_to_name = df_plot.drop_duplicates(subset=["ADseq_"]).set_index("ADseq_")["Name_"].to_dict()
        
        for col_idx, rep in enumerate(reps):
            sub = df_plot[df_plot["rep_"] == rep]
            
            # Top row: z-scored activity
            ax_top = axes[0, col_idx]
            for adseq in adseqs:
                sub_ad = sub[sub["ADseq_"] == adseq]
                ax_top.plot(sub_ad["time_"], sub_ad["z-scored_activity_mean"], 
                            color=adseq_color[adseq], alpha=0.8, marker = 'o', markersize = 5)
                ax_top.fill_between(sub_ad["time_"], 
                                    sub_ad["z-scored_activity_mean"] - sub_ad["z-scored_activity_std"],
                                    sub_ad["z-scored_activity_mean"] + sub_ad["z-scored_activity_std"],
                                    color=adseq_color[adseq], alpha=0.15)
            ax_top.set_ylabel("")
            
            # Bottom row: shifted activity
            ax_bot = axes[1, col_idx]
            for adseq in adseqs:
                sub_ad = sub[sub["ADseq_"] == adseq]
                ax_bot.plot(sub_ad["time_"], sub_ad["shifted_activity_mean"], 
                            color=adseq_color[adseq], alpha=0.8, marker = 'o', markersize = 5)
                ax_bot.fill_between(sub_ad["time_"], 
                                    sub_ad["shifted_activity_mean"] - sub_ad["shifted_activity_std"],
                                    sub_ad["shifted_activity_mean"] + sub_ad["shifted_activity_std"],
                                    color=adseq_color[adseq], alpha=0.15)
            ax_bot.set_xlabel("Time")
            ax_bot.set_ylabel("")
        
        axes_flat = axes.flatten()
        axes_flat[0].set_ylabel("Z-scored Activity")
        axes_flat[3].set_ylabel("Centered Activity")
        axes_flat[0].set_title("Rep 2")
        axes_flat[1].set_title("Rep 3")
        axes_flat[2].set_title("Rep 4")
        axes_flat[3].set_xlabel("")
        axes_flat[5].set_xlabel("")
        
        # Optional: single legend for all ADseqs using Name_
        handles = [plt.Line2D([0], [0], color=adseq_color[ad], lw=2) for ad in adseqs]
        labels = [adseq_to_name[ad] for ad in adseqs]  # map ADseq -> Name_
        fig.legend(handles, labels, bbox_to_anchor=(1.3, 1), loc="upper right", fontsize='small', )
        fig.align_ylabels(axes)
        
        sns.despine()
        plt.tight_layout(pad = 0.1)
        fig.suptitle(group_name, y = 1.07, x = 0.53)
        plt.show()

def plot_control_groups(
    all_activities,
    value_col,
    second_value_col=None,
    controls_path="../../../data/chopTFs_controls.csv",
    reps=(2, 3, 4),
):
    # --- Load + merge ---
    chopTF_controls = pd.read_csv(controls_path)
    all_activities = all_activities.rename(columns={"seq": "ADseq"})
    
    df = pd.merge(
        all_activities,
        chopTF_controls[["ADseq", "Name", "group"]],
        on="ADseq"
    )

    # --- Columns to aggregate ---
    cols_to_agg = [value_col]
    if second_value_col:
        cols_to_agg.append(second_value_col)

    grouped = (
        df.groupby(["ADseq", "rep", "time", "Name", "group"])[cols_to_agg]
        .agg(["mean", "std"])
        .reset_index()
    )

    # flatten columns
    grouped.columns = [
        "_".join(col).strip() if isinstance(col, tuple) else col
        for col in grouped.columns
    ]

    # --- Loop over groups ---
    for group_name in grouped["group_"].unique():
        df_plot = grouped[grouped["group_"] == group_name].sort_values("time_")

        n_rows = 2 if second_value_col else 1
        n_cols = len(reps)

        fig, axes = plt.subplots(
            n_rows, n_cols,
            figsize=(4 * n_cols, 4 * n_rows),
            sharex=True,
            sharey='row' if second_value_col else True,
            dpi=300
        )

        if n_rows == 1:
            axes = axes.reshape(1, -1)

        # --- Colors ---
        adseqs = df_plot["ADseq_"].unique()
        palette = sns.color_palette("tab10", n_colors=len(adseqs))
        adseq_color = dict(zip(adseqs, palette))

        adseq_to_name = (
            df_plot.drop_duplicates("ADseq_")
            .set_index("ADseq_")["Name_"]
            .to_dict()
        )

        # --- Plotting helper ---
        def plot_metric(ax, sub, metric):
            for adseq in adseqs:
                sub_ad = sub[sub["ADseq_"] == adseq]

                mean_col = f"{metric}_mean"
                std_col = f"{metric}_std"

                ax.plot(
                    sub_ad["time_"],
                    sub_ad[mean_col],
                    color=adseq_color[adseq],
                    alpha=0.8
                )

                ax.fill_between(
                    sub_ad["time_"],
                    sub_ad[mean_col] - sub_ad[std_col],
                    sub_ad[mean_col] + sub_ad[std_col],
                    color=adseq_color[adseq],
                    alpha=0.15
                )

        # --- Loop over reps ---
        for col_idx, rep in enumerate(reps):
            sub = df_plot[df_plot["rep_"] == rep]

            # Row 1
            plot_metric(axes[0, col_idx], sub, value_col)
            axes[0, col_idx].set_title(f"Rep {rep}")

            # Row 2 (optional)
            if second_value_col:
                plot_metric(axes[1, col_idx], sub, second_value_col)

        # --- Labels ---
        axes[0, 0].set_ylabel(value_col)
        if second_value_col:
            axes[1, 0].set_ylabel(second_value_col)

        for ax in axes[-1]:
            ax.set_xlabel("Time")

        # --- Legend ---
        handles = [
            plt.Line2D([0], [0], color=adseq_color[ad], lw=2)
            for ad in adseqs
        ]
        labels = [adseq_to_name[ad] for ad in adseqs]

        fig.legend(
            handles, labels,
            bbox_to_anchor=(1.3, 1),
            loc="upper right",
            fontsize="small"
        )

        sns.despine()
        fig.align_ylabels(axes)

        fig.suptitle(group_name, y=1.05)
        plt.tight_layout(pad=0.2)
        plt.show()