import configparser

def read_config(filename, section):
    config = configparser.ConfigParser()
    config.read(filename)
    
    parameters = {}
    
    if section in config:
        for key, value in config[section].items():
            parameters[key] = value
    
    return parameters


# def read_fasta(file_path):
#     """
#     Read a FASTA file and return the sequence as a string.
    
#     Parameters:
#     - file_path (str): Path to the FASTA file.

#     Returns:
#     - str: Sequence data.
#     """
#     sequence = ""
#     with open(file_path, 'r') as file:
#         for line in file:
#             # Skip header lines starting with '>'
#             if not line.startswith('>'):
#                 # Remove newline characters and concatenate the sequence
#                 sequence += line.strip()
#     return sequence


def read_fasta(file_path):
    """
    Read a FASTA file and return the sequence as a string.
    Parameters:
    - file_path (str): Path to the FASTA file.
    Returns:
    - str: Sequence data.
    """
    sequence = ""
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if not line.startswith('>'):
                sequence += line
    return sequence

# def read_fasta(file_path):
#     """
#     Read a FASTA file and return the sequence as a string.

#     Parameters:
#     - file_path (str): Path to the FASTA file.

#     Returns:
#     - dict: Dictionary with sequence data.
#     """
#     sequences = {}
#     current_sequence = ""
#     current_key = None

#     with open(file_path, 'r') as file:
#         for line in file:
#             line = line.strip()
#             if line.startswith('>'):
#                 if current_key is not None:
#                     sequences[current_key] = current_sequence
#                 current_key = line[1:]
#                 current_sequence = ""
#             else:
#                 current_sequence += line

#         if current_key is not None:
#             sequences[current_key] = current_sequence
#     return sequences
    
    
# def read_fasta(file_path):
#     sequences = {}
#     current_description = ""
#     current_sequence = ""

#     with open(file_path, 'r') as fasta_file:
#         for line in fasta_file:
#             line = line.strip()
#             if line.startswith('>'):
#                 # Save the previous sequence if any
#                 if current_description and current_sequence:
#                     sequences[current_description] = current_sequence
#                 # Start a new sequence
#                 current_description = line[1:].strip()  # Remove the ">" and leading/trailing spaces
#                 current_sequence = ""
#             else:
#                 current_sequence += line

#         # Save the last sequence
#         if current_description and current_sequence:
#             sequences[current_description] = current_sequence

#     return sequences

# Example usage:
# fasta_file_path = "..\\output\\paternal_genome.fasta"
# sequences = read_fasta(fasta_file_path)

# Print the result
# for description, sequence in sequences.items():
#     print(f"Description: {description}")
#     print(f"Sequence Length: {len(sequence)}")



def write_fasta(file_path, sequence, description=""):
    with open(file_path, "w") as fasta_file:
        fasta_file.write(f">{description}\n")
        # Write the entire sequence in a single line
        fasta_file.write(sequence + "\n")
    
import time

def process_chunk(function, chunk, result_list):
    start_time = time.time()
    function(**chunk)
    end_time = time.time()
    result_list.append(end_time - start_time)









