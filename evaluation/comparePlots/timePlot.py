import os
import glob
import matplotlib.pyplot as plt
import pandas as pd

def read_elapsed_time(file_path):
    elapsed_times = []

    with open(file_path, 'r') as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith("Elapsed time for simulation"):
            parts = line.split(":")
            elapsed_time = float(parts[1].strip().split()[0])
            elapsed_times.append(elapsed_time / 60.0)  # Convert to minutes

    return elapsed_times

def main(folder_path):
    all_elapsed_times = []

    for file_path in glob.glob(os.path.join(folder_path, '*.log')):
        elapsed_times = read_elapsed_time(file_path)
        all_elapsed_times.append(elapsed_times)

    df = pd.DataFrame(all_elapsed_times)

    # Create a new column for simulation number
    df['Simulation'] = (df.index % 10) + 1

    # Create a new column for condition (300k, 400k, 500k)
    df['Condition'] = df.index // 10
    df['Condition'] = df['Condition'].map({0: 300, 1: 400, 2: 500})

    # Melt the DataFrame to combine 'Simulation' and 'Condition' into one column
    melted_df = pd.melt(df, id_vars=['Simulation', 'Condition'], value_name='Elapsed Time')

    # Create a box plot
    plt.figure(figsize=(12, 6))
    melted_df.boxplot(column='Elapsed Time', by=['Condition', 'Simulation'], showfliers=False)
    
    plt.title('PTA computational time')
    plt.xlabel('Condition')
    plt.ylabel('Elapsed Time (minutes)')
    plt.suptitle('')
    plt.xticks(range(1, 4), ['300k', '400k', '500k'])  # Set custom x-axis labels
    plt.show()

if __name__ == "__main__":
    folder_path = '/path/to/your/folder'  # Replace with the actual path to your folder
    main(folder_path)




folder_path = '..\\timeAnalysis\\PTA\\'
