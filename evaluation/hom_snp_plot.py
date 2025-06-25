import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

csv_file = "/home/bahonar/bahonar_work/SRP178894/haplotypecaller/homo_data.csv"
dp_homo_data = pd.read_csv(csv_file, sep='\t')
homo_data = pd.DataFrame(dp_homo_data)

plt.figure(figsize=(8, 6))
sns.boxplot(x='type', y='cor', hue='range', data=homo_data, palette="Set1", width=0.3, whis=1.5, showfliers=False) 
plt.xlabel('Correlations across repeat simulations')
plt.xticks(rotation=45)
plt.ylabel('Correlations')
plt.title('Correlations of Depth of coverage based on physical distances')
plt.savefig("../test/homo.png", format='png')  # Save the plot
plt.show()