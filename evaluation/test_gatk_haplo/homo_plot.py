import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

pta_rates = [56.36,55.87,54.15,55.37,55.93,55.04,54.35,56.27,56.27,54.96]
mda_rates = [14.88, 19.59,17.34,33.12,39.39,38.13,47.62,35.95,24.82]

df = pd.DataFrame({
    "Method": ["PTA"] * len(pta_rates) + ["MDA"] * len(mda_rates),
    "MisclassificationRate": pta_rates + mda_rates
})

plt.figure(figsize=(8, 6))
ax = sns.boxplot(data=df, x="Method", y="MisclassificationRate", palette="Set2")
sns.stripplot(data=df, x="Method", y="MisclassificationRate", color='black', alpha=0.5, jitter=True)

for method in df['Method'].unique():
    mean_val = df[df['Method'] == method]['MisclassificationRate'].mean()
    x = list(df['Method'].unique()).index(method)
    ax.text(x, mean_val + 1, f"Mean: {mean_val:.1f}%", ha='center', fontsize=10, color='black')

plt.title("Homozygous SNPs Miscalled as Heterozygous")
plt.ylabel("Misclassification Rate (%)")
plt.xlabel("Amplification Method")
plt.tight_layout()
plt.show()