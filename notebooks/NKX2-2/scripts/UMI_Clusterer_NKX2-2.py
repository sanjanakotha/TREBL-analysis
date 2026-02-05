from umi_tools import UMIClusterer
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

directional_clusterer = UMIClusterer(cluster_method="directional")
adjacency_clusterer = UMIClusterer(cluster_method="adjacency")

# Grouped NKX2-2 reads
grouped = pd.read_csv("../../data/NKX2-2_grouped.csv", index_col = 0)

# Need quality ADs since all must be the same length
grouped_AD_qual = grouped[grouped["AD_qual"]]
counts_per_AD = grouped_AD_qual[["AD", "count", "Designed"]].groupby("AD").sum()

# Convert to dictionary
counts_per_AD_dict = counts_per_AD.to_dict()['count']
# Convert all keys to bytes safely (if not already bytes)
counts_per_AD_dict_bytes = {
    (k.encode('utf-8') if isinstance(k, str) else k): v
    for k, v in counts_per_AD_dict.items()
}

# print("Starting directional clustering...")
# directional_clustered_ADs = directional_clusterer(counts_per_AD_dict_bytes, threshold=1)
# directional_clustered_ADs_df = pd.DataFrame(directional_clustered_ADs)
# directional_clustered_ADs_df['group_size'] = directional_clustered_ADs_df.notna().sum(axis=1)
# # Convert all bytes to strings safely
# directional_clustered_ADs_df = directional_clustered_ADs_df.applymap(
#     lambda x: x.decode() if isinstance(x, bytes) else x
# )
# directional_clustered_ADs_df.to_csv("../../data/CC_NX22-2_ADs_directional_clustered.csv")
# print("Done with directional clustering of ADs.")

print("Starting adjacency clustering...")
adjacency_clustered_ADs = adjacency_clusterer(counts_per_AD_dict_bytes, threshold=2)
adjacency_clustered_ADs_df = pd.DataFrame(adjacency_clustered_ADs)
adjacency_clustered_ADs_df['group_size'] = adjacency_clustered_ADs_df.notna().sum(axis=1)
# Convert all bytes to strings safely
adjacency_clustered_ADs_df = adjacency_clustered_ADs_df.applymap(
    lambda x: x.decode() if isinstance(x, bytes) else x
)
adjacency_clustered_ADs_df.to_csv("../../data/CC_NX22-2_ADs_adjacency_clustered_dist_2.csv")
print("Done with adjacency clustering of ADs.")

