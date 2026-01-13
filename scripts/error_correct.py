import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from scripts.preprocess import time_it
from pathlib import Path
import pandas as pd
import matplotlib
import os

matplotlib.use("Agg")  # non-GUI backend for headless plotting

import matplotlib.pyplot as plt
from pathlib import Path
import seaborn as sns
from Bio.Seq import Seq

conda_env='/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/umi_tools/'

def rc_and_swap(preceder: str, post: str, length: int) -> (str, str):
    """
    Given original `preceder` and `post`, return new preceder_rc and post_rc
    for reverse‑complement mode: preceder_rc = RC(post), post_rc = RC(preceder).
    """
    preceder_rc = str( Seq(post).reverse_complement() )
    post_rc     = str( Seq(preceder).reverse_complement() )
    return preceder_rc, post_rc
    
@time_it
def concatenate_fastqs(input_fastq_list, output_dir):
    """
    Concatenate multiple FASTQ files into a single temporary FASTQ in output_dir.
    
    Returns the path to the concatenated file.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    tmp_fastq = output_dir / "tmp_combined.fastq"
    
    # If only one file, just return it
    if isinstance(input_fastq_list, str):
        return input_fastq_list
    
    input_fastq_list = [str(Path(f).resolve()) for f in input_fastq_list]
    with open(tmp_fastq, "wb") as outfile:
        for f in input_fastq_list:
            with open(f, "rb") as infile:
                outfile.write(infile.read())
    return str(tmp_fastq)
    
def run_whitelist_on_concat_domains(fastq_path, output_dir, set_cell_number=None, prefix=None):
    """
    Run umi_tools whitelist on a simplified FASTQ (concatenated barcodes only).

    Args:
        fastq_path (str or Path): Path to the FASTQ file containing extracted barcodes.
        output_dir (str or Path): Directory to save whitelist outputs.
        prefix (str, optional): Prefix for naming log and output files (default = FASTQ basename).
    """
    fastq_path = Path(fastq_path).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(exist_ok=True)

    if prefix is None:
        prefix = fastq_path.stem  # e.g. "HAR_AD_filtered_barcode"

    log_path = output_dir / f"{prefix}_whitelist.log"
    whitelist_path = output_dir / f"{prefix}_whitelist.txt"
    plot_prefix = output_dir / f"{prefix}_plots"

    cmd = [
        "conda", "run", "-p", conda_env,
        "umi_tools", "whitelist",
        "--knee-method=density",
        "--stdin", str(fastq_path),
        "--bc-pattern", "(?P<umi_1>.{1})(?P<cell_1>.*)",
        "--extract-method=regex",
        "--method=reads",
        "--log", str(log_path),
        "--stdout", str(whitelist_path),
        "--plot-prefix", str(plot_prefix)
    ]

    if set_cell_number is not None:
        print("Using custom cell number.")
        cmd.extend(["--set-cell-number", str(set_cell_number)])

    print(f"Running umi_tools whitelist on {fastq_path.name} ...")
    subprocess.run(cmd, check=True)
    
    print(f"Whitelist complete.\n- Log: {log_path}\n- Output: {whitelist_path}\n- Plots: {plot_prefix}_*.png")

    return {
        "log": log_path,
        "whitelist": whitelist_path,
        "plot_prefix": plot_prefix
    }

def convert_txt_to_whitelist_mapping_df_from_path(whitelist_path):
    """
    Read a umi_tools whitelist text file (single concatenated FASTQ) and
    build a mapping DataFrame identical to the old per-barcode function.

    Parameters:
        whitelist_path (str or Path): Path to the whitelist txt file.
        reverse_complement (bool): Whether to reverse complement sequences.

    Returns:
        pd.DataFrame: Mapping DataFrame with columns 'original' and 'canonical'.
    """
    

    whitelist_path = Path(whitelist_path)
    
    # Assume umi_tools output has canonical, collapsed, largest_count, counts
    wl_df = pd.read_csv(
        whitelist_path, sep="\t", header=None,
        names=["canonical", "collapsed", "largest_count", "counts"]
    )

    def maybe_revcomp(seq):
        return seq
        
    all_rows = []
    for _, row in wl_df.iterrows():
        canonical = maybe_revcomp(row["canonical"])
        all_rows.append({"original": canonical, "canonical": canonical})
        
        if pd.notna(row["collapsed"]):
            for c in row["collapsed"].split(","):
                all_rows.append({"original": maybe_revcomp(c), "canonical": canonical})

    mapping_df = pd.DataFrame(all_rows)
    return mapping_df


# @time_it
# def run_single_whitelist(barcode, input_fastq, output_dir, reverse_complement):
#     print(f"Generating whitelist for {barcode.name} (RC={reverse_complement})…")

#     if reverse_complement:
#         preceder_rc, post_rc = rc_and_swap(barcode.preceder, barcode.post, barcode.length)
#         bc_pattern = f".*{preceder_rc}(?P<cell_1>.{{{barcode.length}}}){post_rc}(?P<umi_1>.{{1}}).*"
#     else:
#         bc_pattern = f".*{barcode.preceder}(?P<cell_1>.{{{barcode.length}}}){barcode.post}(?P<umi_1>.{{1}}).*"

#     log_path       = f"{output_dir}/{barcode.name}_extract.log"
#     whitelist_path = f"{output_dir}/{barcode.name}_whitelist.txt"

#     cmd = [
#         "conda", "run", "-p", conda_env, "umi_tools", "whitelist",
#         "--stdin", input_fastq,
#         "--bc-pattern", bc_pattern,
#         "--extract-method", "regex",
#         "--method", "reads",
#         "--log", log_path,
#         "--stdout", whitelist_path
#     ]

#     subprocess.run(cmd, check=True)
#     return f"Whitelist generated for {barcode.name} at {whitelist_path}"
    

# def run_whitelist_parallel(barcodes, input_fastq, output_dir, reverse_complement, max_workers=8):
#     """
#     Run umi_tools whitelist in parallel for multiple Barcode objects.

#     Parameters:
#         barcodes (list[Barcode]): List of Barcode objects.
#         input_fastq (str or list[str]): Path(s) to input FASTQ files.
#         output_dir (str): Directory to save logs and whitelists.
#         max_workers (int): Number of parallel processes.
#     """
#     # Concatenate FASTQs if it's a list
#     input_fastq_path = concatenate_fastqs(input_fastq, output_dir)

#     results = []
#     with ThreadPoolExecutor(max_workers=max_workers) as executor:
#         future_to_barcode = {executor.submit(run_single_whitelist, bc, input_fastq_path, output_dir, reverse_complement): bc for bc in barcodes}
#         for future in as_completed(future_to_barcode):
#             try:
#                 result = future.result()
#                 results.append(result)
#                 print(result)
#             except subprocess.CalledProcessError as e:
#                 bc = future_to_barcode[future]
#                 print(f"Error processing {bc.name}: {e}")
#     return results

# def convert_txt_to_whitelist_mapping_df(barcode, output_dir, reverse_complement):
#     """
#     Read the whitelist file for a single Barcode object and build a mapping
#     DataFrame from collapsed barcodes to canonical barcodes.
#     Reverse complement barcodes if reverse_complement=True.

#     Parameters:
#         barcode (Barcode): The Barcode object, must have attribute 'name'.
#         output_dir (str or Path): Directory where the whitelist file is saved.
#         reverse_complement (bool): Whether to reverse complement barcodes.

#     Returns:
#         pd.DataFrame: Mapping DataFrame with columns 'original' and 'canonical'.
#     """
#     output_dir = Path(output_dir)
#     whitelist_path = output_dir / f"{barcode.name}_whitelist.txt"
    
#     wl_df = pd.read_csv(
#         whitelist_path, sep="\t", header=None,
#         names=["canonical", "collapsed", "largest_count", "counts"]
#     )

#     def maybe_revcomp(seq):
#         if pd.isna(seq):
#             return seq
#         return str(Seq(seq).reverse_complement()) if reverse_complement else seq

#     all_rows = []
#     for _, row in wl_df.iterrows():
#         canonical = maybe_revcomp(row["canonical"])
#         all_rows.append({"original": canonical, "canonical": canonical})
        
#         if pd.notna(row["collapsed"]):
#             for c in row["collapsed"].split(","):
#                 all_rows.append({"original": maybe_revcomp(c), "canonical": canonical})
    
#     mapping_df = pd.DataFrame(all_rows)
#     return mapping_df

def summarize_whitelist(barcode, output_dir, reverse_complement, fig=None, axs=None, figsize=(16, 4), dpi=300):
    """
    Summarize a umi_tools whitelist for a single Barcode and plot 4 summary plots side by side:
        1. Bar chart: original sequences vs canonical sequences
        2. Histogram of group sizes
        3. Scatter: largest member vs group size
        4. Scatter: largest member vs sum of smaller members

    Parameters:
        barcode (Barcode): Barcode object whose whitelist to summarize.
        output_dir (str or Path): Directory containing the whitelist output.
        fig (matplotlib.figure.Figure, optional): Pre-created figure. If None, new figure is created.
        axs (array-like, optional): Pre-created matplotlib axes to plot on. If None, new axes are created.
        figsize (tuple): Figure size if creating a new figure.
        dpi (int): Figure DPI if creating a new figure.

    Returns:
        matplotlib.figure.Figure: Figure object used for plotting
        array-like: Axes array used for plotting
    """
    
    output_dir = Path(output_dir)
    whitelist_path = output_dir / f"{barcode.name}_whitelist.txt"
    
    wl_df = pd.read_csv(
        whitelist_path, sep="\t", header=None,
        names=["canonical", "collapsed", "largest_count", "counts"]
    )   
    
    # Number of sequences per canonical
    wl_df["collapsed_list"] = wl_df["collapsed"].fillna("").apply(lambda x: x.split(",") if x else [])
    wl_df["num_merged"] = wl_df["collapsed_list"].apply(len)
    
    def parse_counts(row):
        if pd.isna(row["counts"]) or str(row["counts"]).strip() == "":
            return row["largest_count"], 0
        rest = sum(int(c) for c in str(row["counts"]).split(",") if c.strip())
        largest = row["largest_count"] - rest
        return largest, rest
    
    wl_df[["largest_count", "rest_count"]] = wl_df.apply(parse_counts, axis=1, result_type="expand")
    
    # Total counts for bar chart
    total_original_seqs = wl_df["num_merged"].sum() + len(wl_df)
    total_canonical_seqs = len(wl_df)

    # Create figure/axes if not provided
    if fig is None or axs is None:
        fig, axs = plt.subplots(1, 4, figsize=figsize, dpi=dpi)

    # Panel 1: Bar chart of original vs canonical sequences
    axs[0].bar(["Before", "After"], [total_original_seqs, total_canonical_seqs],
               color=[sns.color_palette('Paired')[0], sns.color_palette('Paired')[1]])
    axs[0].set_ylabel("Number of sequences")
    # Annotate counts above each bar
    for x, y in zip(["Before", "After"], [total_original_seqs, total_canonical_seqs]):
        axs[0].text(x, y + max([total_original_seqs, total_canonical_seqs])*0.02, f"{y:,}", 
                    ha='center', va='bottom', fontsize='medium', weight='bold')

    # Panel 2: Histogram of group sizes
    axs[1].hist(wl_df["num_merged"] + 1, bins=30, edgecolor='black')
    axs[1].set_xlabel("Group size")
    axs[1].set_ylabel("Frequency")
    
    # Panel 3: Scatter - largest member vs group size
    axs[2].scatter(wl_df["largest_count"], wl_df["num_merged"] + 1, alpha=0.6)
    axs[2].set_xlabel("Reads of largest")
    axs[2].set_ylabel("Group size")
    
   # Panel 4: Scatter - largest member vs sum of smaller members
    axs[3].scatter(wl_df["largest_count"], wl_df["rest_count"], alpha=0.6)
    axs[3].set_xlabel("Reads of largest")
    axs[3].set_ylabel("Summed merged reads")
    axs[3].set_xscale('log')
    axs[3].set_yscale('log')
    
    # Save current axis limits
    xlims = axs[3].get_xlim()
    ylims = axs[3].get_ylim()
    
    # Determine line range within current view
    line_min = min(xlims[0], ylims[0])
    line_max = max(xlims[1], ylims[1])
    
    # Plot y=x line
    axs[3].plot([line_min, line_max], [line_min, line_max], 'r--', alpha=0.7)
    
    # Restore axis limits
    axs[3].set_xlim(xlims)
    # axs[3].set_ylim(ylims)
        
    fig.tight_layout(pad=2)
    
    return fig, axs

def plot_all_whitelists(barcodes, output_dir, reverse_complement, n_cols=4, dpi=300):
    """
    Plot summaries of multiple barcodes in a grid layout using summarize_whitelist.

    Parameters
    ----------
    barcodes : list[Barcode]
        List of Barcode objects.
    output_dir : str
        Directory where whitelist outputs are saved.
    n_cols : int, optional
        Number of panels per row. Default is 4.
    dpi : int, optional
        Figure DPI. Default is 300.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object containing all barcode summaries.
    """
    n_barcodes = len(barcodes)

    # Create figure with n_barcodes rows, n_cols columns
    fig, axs = plt.subplots(n_barcodes, n_cols, figsize=(4 * n_cols, 4 * n_barcodes), dpi=dpi)

    # Ensure axs is 2D array for consistent indexing
    if n_barcodes == 1:
        axs = axs.reshape(1, -1)

    # Loop over barcodes
    for i, bc in enumerate(barcodes):
        # Pass the row of axes to summarize_whitelist
        fig, _ = summarize_whitelist(
            barcode=bc,
            output_dir=output_dir,
            fig=fig,
            axs=axs[i],
            reverse_complement = reverse_complement
        )

        # Add row-level label for the barcode
        fig.text(
            -0.02,  # x-position (left margin)
            (n_barcodes - i - 0.5) / n_barcodes,  # y-position (center of row)
            f"{bc.name}",
            fontsize='large',
            weight='bold',
            rotation=90,
            va='center'
        )

    sns.despine()
    plt.tight_layout(pad=2)

    return fig
