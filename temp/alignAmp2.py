import re

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

# Function to count paternal and maternal amplicons that contain the specific point
def count_amplicons_at_point(amplicons, genomic_point):
    pat_count = 0
    mat_count = 0

    # Iterate through the amplicons to check if they cover the specified point
    for name, start, end, _ in amplicons:
        if start <= genomic_point <= end:
            if 'patSeq' in name:
                pat_count += 1
            elif 'matSeq' in name:
                mat_count += 1

    return pat_count, mat_count

# Function to calculate the ratio pat/(pat + mat)
def calculate_pat_mat_ratio(pat_count, mat_count):
    if pat_count + mat_count == 0:
        return None  # No amplicons cover the point
    return pat_count / (pat_count + mat_count)

# Filepath to your FASTA file
fasta_file = '../../test/subset.fasta'  # Replace with your actual file path

# Define the genomic point of interest
genomic_point = 75000  # Replace with the point you are interested in

# Extract amplicons from the FASTA file
amplicons = extract_amplicons_from_fasta(fasta_file)

# Count the number of paternal and maternal amplicons that contain the specific point
pat_count, mat_count = count_amplicons_at_point(amplicons, genomic_point)

# Calculate the ratio pat/(pat + mat)
ratio = calculate_pat_mat_ratio(pat_count, mat_count)

# Print the results
print(f"At genomic point {genomic_point}:")
print(f"Paternal amplicons: {pat_count}")
print(f"Maternal amplicons: {mat_count}")
if ratio is not None:
    print(f"Ratio (pat / (pat + mat)): {ratio:.2f}")
else:
    print("No amplicons cover the specified point.")
