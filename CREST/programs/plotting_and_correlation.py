import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr

df = pd.read_csv("CREST/results/tables/007_JAK1_rotatable_bonds.csv")

metrics = ["rotatable_bonds",
    "num_rings",
    "num_aromatic_rings",
    "num_heavy_atoms",
    "mol_weight",
    "num_HBD",
    "num_HBA"
]
correlation_results = []
for metric in metrics:
    data_clean = df[[metric, "Sconf"]].dropna()


    plt.figure(figsize=(6, 4))
    sns.regplot(data=data_clean, x=metric, y="Sconf", ci=None)
    plt.title(f"Sconf vs. {metric}")
    plt.xlabel(metric)
    plt.ylabel("Sconf")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"CREST/results/plots/Sconf_vs_{metric}.png")
    plt.show()

    pearson_r, pearson_p = pearsonr(data_clean[metric], data_clean["Sconf"])
    spearman_r, spearman_p = spearmanr(data_clean[metric], data_clean["Sconf"])


    print(f"{metric}")
    print(f"  Pearson  r = {pearson_r:.3f} (p = {pearson_p:.3g})")
    print(f"  Spearman r = {spearman_r:.3f} (p = {spearman_p:.3g})")
    print("-" * 40)
    correlation_results.append({
        "metric": metric,
        "pearson_r": pearson_r,
        "pearson_p": pearson_p,
        "spearman_r": spearman_r,
        "spearman_p": spearman_p
    })


corr_df = pd.DataFrame(correlation_results)
corr_df.to_csv("CREST/results/tables/Sconf_correlation_stats.csv", index=False)