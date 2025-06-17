import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import FastaWriter
import os
import pickle
#sys.path.append('/home/bahonar/simulation/SingleCellSim')

####################################################################################

def GenerateNewAmp(refSeq, A, A_i, t, delta_t, main_dtype, lMin, lMax, parent, Gamma, beta, Theta, exclude, n_i):
    A_c = np.zeros(n_i, dtype=main_dtype)
    refSeq_length = len(refSeq)
    ## make new amplicons
    # Specify the direction of A_c
    A_c['direction'] = -A_i['direction']
    A_c['parent'] = parent
    A_c['source'] = A_i['source']
    # Specify start time, start position, and maximum length of A_c
    A_c['startTime'] = np.random.uniform(t, t + delta_t, n_i)
    A_c['maxLength'] = np.random.randint(lMin, lMax, n_i)
    # if exclude:
        # accessibleRegions = computeAccessibleRegions(A, Gamma, parent)
        # A_c['startPos'] = np.random.choice(accessibleRegions, n_i)
    # else:
    A_c['startPos'] = np.random.choice(np.arange(A_i['startPos'], A_i['endPos'], A_i['direction']), n_i)
    end = A_c['startPos'] + ((A_c['maxLength']+1) * A_c['direction'])
    valid_indices = np.where(
        (0 <= A_c['startPos']) & 
        (A_c['startPos'] < refSeq_length) & 
        (0 <= end) & 
        (end < refSeq_length) &
        (abs(end - A_c['startPos']) <= abs(A_i['startPos'] - A_c['startPos']))
    )
    #valid_indices = np.where((0 <= A_c['startPos']) & (A_c['startPos'] < refSeq_length) & (0 <= end) & (end < refSeq_length))
    A_c = A_c[valid_indices]
    n_i = len(valid_indices[0])
    end = A_c['startPos'] + (A_c['maxLength'] * A_c['direction'])
    if len(A_c) != 0:
        updateErrors(A,A_c,parent,end,beta,n_i,refSeq)
        extendAmplicon(A_c, t, delta_t, Theta)
        return A_c
    else:
        return np.array([], dtype=main_dtype)

#####################################################################################

def extendAmplicon(A_i, t, delta_t, Theta):
    A_i['endPos'] = (A_i['startPos'] + Theta * (t + delta_t - A_i['startTime']) * A_i['direction']).astype(int)
    mask = ((A_i['endPos'] - A_i['startPos']) * A_i['direction']) >= A_i['maxLength']
    A_i['endPos'][mask] = (A_i['startPos'] + ((A_i['maxLength'] +1)* A_i['direction']))[mask]
    A_i['released'][mask] = True
    
    return A_i

####################################################################################

# def computeAccessibleRegions(A, Gamma, parent):
#     excludedRegions = np.zeros_like(A, dtype=bool)

#     for A_c in A:
#         mask = (A_c['parent'] == parent) & ~A_c['released']
#         start = A_c['endPos'][mask]
#         end = (A_c['endPos'] - (Gamma * A_c['direction']))[mask]
#         excludedRegions[mask] = (start > end)

#     # Accessible regions are the difference between the entire region and excluded regions
#     accessibleRegions = np.where(~excludedRegions)
    
#     return accessibleRegions
# ####################################################################################
def getAlt(currentNucleotide):
    possibleNucleotides = [nt for nt in 'acgtACGT' if nt.lower() != currentNucleotide.lower()]    
    # Ensure the case of the alternate nucleotide matches the case of the current nucleotide
    alternate_nucleotide = np.random.choice(possibleNucleotides).upper() if currentNucleotide.isupper() else np.random.choice(possibleNucleotides).lower()    
    return alternate_nucleotide

####################################################################################

def updateErrors(A, A_c, parent, end, beta, n_i, refSeq):
    # 1. Inherit errors from parent
    parent_errors = A[parent]['errors']

    # Initialize inherited error lists for each amplicon
    inherited_errors = [[] for _ in range(n_i)]

    if len(parent_errors) > 0:
        for i, (start_pos, end_pos, direction) in enumerate(zip(A_c['startPos'], end, A_c['direction'])):
            errors_in_range = []
            for pos, nt in parent_errors:
                if direction == -1 and end_pos <= pos <= start_pos:
                    errors_in_range.append((pos, nt))
                elif direction == 1 and start_pos <= pos <= end_pos:
                    errors_in_range.append((pos, nt))
            inherited_errors[i] = errors_in_range

    # 2. Add new errors based on binomial
    error_positions_count = np.random.binomial(A_c['maxLength'], beta, n_i)
    temp_start = np.minimum(A_c['startPos'], end)
    temp_end = np.maximum(A_c['startPos'], end)

    for i in range(n_i):
        new_errors = []
        if error_positions_count[i] > 0:
            positions = np.random.randint(temp_start[i], temp_end[i] + 1, size=error_positions_count[i])
            new_errors = [(pos, getAlt(refSeq[pos])) for pos in positions]

        # Append new errors to inherited
        inherited_errors[i].extend(new_errors)

    A_c['errors'] = inherited_errors
    return A_c

####################################################################################
# MDASimulation: this function get the positions of amplicons that is produced in the amplification simulation
# Theta: the number of nucleotides that is made by polymerase in MDA time unit based on polymerase speed
# Gamma: the number of nucleotides that is displaced 
def MDASimulation(patSeq, matSeq, Theta=12000, Gamma=50, DNACoef=400,
                            lMin=2000, lMax=70000, Lambda=0.0001,
                            delta_t=0.1, beta=0.0001, exclude=False,
                            saveInterval=15, output_folder="output", resume=False,
                            amp_depth = 15, template = False):
    P = len(patSeq)
    M = len(matSeq)
    main_dtype = np.dtype([
        ('startPos', int),
        ('endPos', int),
        ('maxLength', int),
        ('direction', int),
        ('parent', int),
        ('errors', 'O'),
        ('released', bool),
        ('startTime', float),
        ('source', 'U1')
    ])
    # check_DNA = (P + M) * 2
    check_DNA = 0
    initial_DNA = P
    final_DNA = P * amp_depth
    A = np.array([
        (0, M-1, M, +1, -1, [], True, 0.0, 'M'),
        (M-1, 0, M, -1, -1, [], True, 0.0, 'M'),
        (0, P-1, P, +1, -1, [], True, 0.0, 'P'),
        (P-1, 0, P, -1, -1, [], True, 0.0, 'P')
    ], dtype=main_dtype)
    t = 0
    total_maxLength = 0 
    num_amplicons = 0
    At = A
    if template:
        while check_DNA <= final_DNA:
            # Generate new amplicons for all parents using vectorized approach
            n_i_values = np.random.poisson(np.abs(Lambda * delta_t * (((A['endPos'] - A['startPos']) * A['direction'])+1) ))
            # print(n_i_values)
            # print(sum(n_i_values[4:]))
            newAmplicons = np.concatenate([GenerateNewAmp(patSeq if A[parent]['source'] == 'P' else matSeq,
                                                          A, A[parent], t,  delta_t, main_dtype, lMin, lMax, parent, Gamma, beta, Theta, exclude, n_i)
                                           for parent, n_i in enumerate(n_i_values)
            ])
            unreleased_mask = ~A['released']
            A[unreleased_mask] = extendAmplicon(A[unreleased_mask], t, delta_t, Theta)
            # Concatenate the new amplicons to A using vectorized approach
            A = np.concatenate((A, newAmplicons))
            t += delta_t
            total_maxLength += np.sum(newAmplicons['maxLength'])
            # print(total_maxLength)
            num_amplicons += len(newAmplicons)
            # print(t)
            # print('everything well done')
            check_DNA += np.sum(newAmplicons['maxLength']) 
            #print(check_DNA)
    else:
        while check_DNA <= final_DNA:
            # Generate new amplicons for all parents using vectorized approach
            n_i_values = np.random.poisson(np.abs(Lambda * delta_t * (At['endPos'] - At['startPos']) * At['direction']))
            #print(n_i_values)
            newAmplicons = np.concatenate([GenerateNewAmp(patSeq if At[parent]['source'] == 'P' else matSeq,
                                                          At, At[parent], t, delta_t, main_dtype, lMin, lMax, parent, Gamma, beta, Theta, exclude, n_i)
                                           for parent, n_i in enumerate(n_i_values)
            ])
            unreleased_mask = ~A['released']
            A[unreleased_mask] = extendAmplicon(A[unreleased_mask], t, delta_t, Theta)
            # Concatenate the new amplicons to A using vectorized approach
            A = np.concatenate((A, newAmplicons))
            t += delta_t
            total_maxLength += np.sum(newAmplicons['maxLength'])
            #print(total_maxLength)
            num_amplicons += len(newAmplicons)
            #print(t)
            check_DNA += np.sum(newAmplicons['maxLength'])
            #print(check_DNA)
            
    #print(t)
    # print(num_amplicons)
    average_maxLength = total_maxLength / num_amplicons
    # print("Average maxLength:", average_maxLength)
    # subseting amplicons for the desired depth 
    weights  = [entry['maxLength']/total_maxLength for entry in A[np.arange(4,len(A))]]
    subset = int((amp_depth*initial_DNA)/average_maxLength)
    # print(subset)
    
    # print(subset_percentage)
    # num_elements_to_subset = int(len(A) * subset_percentage)
    # print(num_elements_to_subset)
    # print(len(A))
    random_indices = np.random.choice(np.arange(4,len(A)), subset, replace=False)#, p= weights)
    selected_indices = np.concatenate(([0, 1, 2, 3], random_indices))
    subsetA = A[selected_indices]
    # print(len(subsetA))
    total_length_maxLength = np.sum(subsetA['maxLength'])
    print(total_length_maxLength)
    outputFilename = os.path.join(output_folder, "amplicons.pkl")
    with open(outputFilename, 'wb') as f:  
        pickle.dump(A, f)
    return A 

####################################################################################
def applyErrors(sequence, errors, startPos):
    """
    Applies substitution errors to the given sequence.

    Parameters:
    - sequence (str): The original DNA sequence.
    - errors (list of tuples): Each tuple is (position, alt_base), where position is absolute.
    - startPos (int): The absolute start position of the sequence.

    Returns:
    - str: The mutated sequence.
    """
    seq_list = list(sequence)
    seq_len = len(seq_list)

    for pos, alt_base in errors:
        idx = pos - startPos  # Convert absolute to relative position
        if 0 <= idx < seq_len:
            seq_list[idx] = alt_base
        else:
            print(f"Warning: position {pos} (relative {idx}) out of bounds for sequence starting at {startPos} with length {seq_len}")
    
    return ''.join(seq_list)

    
####################################################################################
# subsetAmpliconSaveTofFASTA: this function first subset a desired percentage of amplicons randomly
# then convert them to a fasta file that is ready for the next step, sequencing by Art

def subsetAmpliconSaveToFASTA(amplicons, patSeq, matSeq, output_folder="output"):
    """
    Converts a list of amplicons with errors into a FASTA file.

    Parameters:
    - amplicons (list): List of amplicon dictionaries containing keys:
        'startPos', 'maxLength', 'direction', 'source', 'errors'
    - patSeq (str): The paternal reference sequence.
    - matSeq (str): The maternal reference sequence.
    - output_folder (str): The folder where the FASTA file will be saved.
    """

    records = []

    for i, amplicon in enumerate(amplicons):
        length = amplicon['maxLength']
        direction = amplicon['direction']
        source = amplicon['source']
        errors = amplicon['errors']

        if direction == -1:
            startPos = amplicon['startPos'] + length * direction
            endPos = amplicon['startPos']
        else:
            startPos = amplicon['startPos']
            endPos = amplicon['startPos'] + length * direction

        # Select the correct reference sequence
        refSeq = patSeq if source == 'P' else matSeq
        seqName = "patSeq" if source == 'P' else "matSeq"

        # Extract the subsequence and apply errors
        ampliconSeq = refSeq[startPos:endPos]
        ampliconSeq = applyErrors(ampliconSeq, errors, startPos)

        # Reverse complement if needed
        if direction == -1:
            ampliconSeq = str(Seq(ampliconSeq).reverse_complement())
        else:
            ampliconSeq = str(Seq(ampliconSeq))

        # Create SeqRecord
        record_id = f"{seqName}_{direction}_{i}"
        record_description = f"Start: {startPos+1}, End: {endPos}"
        record = SeqRecord(Seq(ampliconSeq), id=record_id, description=record_description)

        records.append(record)

    # Make sure the output directory exists
    os.makedirs(output_folder, exist_ok=True)
    outputFilename = os.path.join(output_folder, "subset.fasta")

    # Write sequences to FASTA with 60 characters per line
    with open(outputFilename, "w") as fasta_file:
        writer = FastaWriter(fasta_file, wrap=60)
        writer.write_file(records)


#####################################################################################
# def test_applyErrors_and_amplicons():
#     patSeq = "ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG"
#     errors = [(10, 'A'), (15, 'T')]

#     amplicon = {
#         'startPos': 5,
#         'maxLength': 10,
#         'direction': 1,
#         'source': 'P',
#         'errors': errors
#     }

#     start = amplicon['startPos']
#     end = start + amplicon['maxLength']
#     ref = patSeq[start:end]
#     expected = list(ref)
#     if 10 - start < len(expected):
#         expected[10 - start] = 'A'
#     if 15 - start < len(expected):
#         expected[15 - start] = 'T'
#     expected = ''.join(expected)

#     result = applyErrors(ref, amplicon['errors'], start)
#     assert result == expected, f"Expected {expected}, got {result}"

#     print("✅ Test passed.")

# patSeq = "ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG"
# matSeq = "TGACTGACTGACTGACTGACTGACTGACTGACTGACTGAC"

# errors = [(10, 'A'), (15, 'G')]

# amplicons = [
#     {
#         'startPos': 5,
#         'maxLength': 10,
#         'direction': 1,
#         'source': 'P',
#         'errors': errors
#     },
#     {
#         'startPos': 8,
#         'maxLength': 10,
#         'direction': -1,
#         'source': 'M',
#         'errors': errors
#     }
# ]

# def test_single_amplicon(amplicon, patSeq, matSeq):
#     length = amplicon['maxLength']
#     direction = amplicon['direction']

#     if direction == -1:
#         startPos = amplicon['startPos'] + length * direction
#         endPos = amplicon['startPos']
#     else:
#         startPos = amplicon['startPos']
#         endPos = amplicon['startPos'] + length * direction

#     refSeq = patSeq if amplicon['source'] == 'P' else matSeq

#     # Slice and apply errors
#     ampliconSeq = refSeq[startPos:endPos]
#     ampliconSeq = applyErrors(ampliconSeq, amplicon['errors'], startPos)

#     # Reverse complement if needed
#     if direction == -1:
#         ampliconSeq = str(Seq(ampliconSeq).reverse_complement())

#     print(f"🧬 Final amplicon sequence:\n{ampliconSeq}")

# test_single_amplicon(amplicons[0], patSeq, matSeq)
