
import re
import matplotlib.pyplot as plt

# Function to extract amplicons from a FASTA file
def extract_amplicons_from_fasta(fasta_file):
    amplicons = []
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Extract sequence name, start, and end using regular expressions
                match = re.search(r'>(\S+).*Start:\s*(\d+),\s*End:\s*(\d+)', line)
                if match:
                    seq_name = match.group(1)
                    start = int(match.group(2))
                    end = int(match.group(3))
                    length = end - start
                    amplicons.append((seq_name, start, end, length))
    return amplicons

# Function to filter amplicons based on region of interest
def filter_amplicons_by_region(amplicons, region_start, region_end):
    filtered_amplicons = []
    for amplicon in amplicons:
        _, start, end, _ = amplicon
        # Check if the amplicon overlaps with the region of interest
        if start <= region_end and end >= region_start:
            filtered_amplicons.append(amplicon)
    return filtered_amplicons

# Function to plot the amplicons within the region of interest
def plot_amplicons_with_focus(amplicons, region_start, region_end):
    # Sort amplicons by length (longest first)
    amplicons_sorted = sorted(amplicons, key=lambda x: x[3], reverse=True)

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot each amplicon, use different colors for patSeq and matSeq
    for i, (name, start, end, length) in enumerate(amplicons_sorted):
        color = 'blue' if 'matSeq' in name else 'red'  # Assign color based on sequence name
        # Plot the amplicon only within the region of interest
        ax.plot([max(start, region_start), min(end, region_end)], [i, i], color=color, linewidth=2, label=name if i == 0 else "")

    # Set axis limits to zoom in on the region of interest
    ax.set_xlim(region_start, region_end)
    ax.set_ylim(-1, len(amplicons_sorted))
    ax.set_xlabel('Genomic Position')
    ax.set_ylabel('Amplicon')

    # Add title
    plt.title(f"Amplicons containing region {region_start}-{region_end}")
    plt.grid(True)

    # Add a legend to distinguish between matSeq and patSeq
    blue_patch = plt.Line2D([0], [0], color='blue', lw=2, label='matSeq')
    red_patch = plt.Line2D([0], [0], color='red', lw=2, label='patSeq')
    ax.legend(handles=[blue_patch, red_patch])

    # Show the plot
    plt.show()

# Filepath to your FASTA file
fasta_file = '../../test/PTA_subset.fasta'  # Replace with your actual file path

# Define the region of interest
region_start = 50000  # Replace with your desired start position
region_end = 52500   # Replace with your desired end position

# Extract amplicons from the FASTA file
amplicons = extract_amplicons_from_fasta(fasta_file)

# Filter amplicons that overlap with the region of interest
filtered_amplicons = filter_amplicons_by_region(amplicons, region_start, region_end)

# Plot the filtered amplicons, focusing on the region
if filtered_amplicons:
    plot_amplicons_with_focus(filtered_amplicons, region_start, region_end)
else:
    print("No amplicons found in the specified region.")
