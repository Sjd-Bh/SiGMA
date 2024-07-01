import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Read data from the file
data = pd.read_csv("..\\..\\outputs\\output.txt", delimiter='\t')

lowest_pos = data['Pos'].min()
highest_pos = data['Pos'].max()

minimum_distance = 1000000  # You can adjust this value according to your requirement
numberOfSelectedPos = 20


# Generate a random position within the specified range
def generate_random_position():
    return np.random.randint(lowest_pos, highest_pos + 1)

positions = [generate_random_position()]

# Keep adding positions that meet the condition
while len(positions) < numberOfSelectedPos:
    # Generate a random position
    random_pos = generate_random_position()
    
    # Check if the random position has minimum distance 'n' from other selected positions
    if all(abs(random_pos - pos) >= minimum_distance for pos in positions):
        positions.append(random_pos)

positions.sort()
end_positions = np.array(positions) + minimum_distance

vaf_ranges_dict = {}

# Iterate through data and store the 'VAF' column for each range
for i in range(len(positions)):
    start_pos = positions[i]
    end_pos = end_positions[i]
    
    # Find 'VAF' values in the specified range
    vaf_values_in_range = data.loc[(data['Pos'] >= start_pos) & (data['Pos'] <= end_pos), 'VAF']
    
    # Store the 'VAF' values in the dictionary only if the range is not empty
    if not vaf_values_in_range.empty:
        range_name = f"{start_pos},{end_pos}"
        vaf_ranges_dict[range_name] = vaf_values_in_range

# Plot a boxplot for the non-empty ranges
if vaf_ranges_dict:
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=pd.DataFrame(vaf_ranges_dict), width=0.5)
    plt.xlabel("Ranges")
    plt.ylabel("VAF")
    plt.title("Boxplot of VAF in Each Range")
    
    # Rotate x-axis labels for better visibility
    plt.xticks(rotation=45, ha="right")  # You can adjust the rotation angle
    
    plt.show()

    # Calculate mean and variance for each range and save the results
    result_dict = {}
    for range_name, vaf_values in vaf_ranges_dict.items():
        mean_vaf = vaf_values.mean()
        var_vaf = vaf_values.var()
        result_dict[range_name] = {'Mean_VAF': mean_vaf, 'Var_VAF': var_vaf}

    # Convert the result to a DataFrame and save it
    result_df = pd.DataFrame(result_dict).T
    result_df.to_csv("..\\VAFanalysis\\PTA\\stat_SRR843_1Mb.csv", index_label='Range')

else:
    print("No non-empty ranges to plot.")
    
    












    