import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def load_results(file_path):
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            sample = fields[0]
            values = {field.split('=')[0]: float(field.split('=')[1]) for field in fields[1:]}
            values['Sample'] = sample
            data.append(values)
    df = pd.DataFrame(data)
    return df

# Example:
prosolo_df = load_results('/home/bahonar/bahonar_work/SRP178894/benchmark/prosolo_summary.tsv')
sccaller_df = load_results('/home/bahonar/bahonar_work/SRP178894/benchmark/sccaller_summary.tsv')

# Check
print(prosolo_df.head())




# Add a 'Tool' column to distinguish
prosolo_df['Tool'] = 'ProSolo'
sccaller_df['Tool'] = 'SCcaller'

# Combine both
combined_df = pd.concat([prosolo_df, sccaller_df], ignore_index=True)

plt.figure(figsize=(8,6))
sns.boxplot(data=combined_df, x='Tool', y='F1')
plt.title('F1-Score Comparison')
plt.ylabel('F1-Score')
plt.xlabel('Variant Caller')
plt.grid(True)
plt.show()


plt.figure(figsize=(8,6))
sns.scatterplot(data=combined_df, x='Recall', y='Precision', hue='Tool', style='Tool', s=100)
plt.title('Precision vs Recall')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.grid(True)
plt.show()