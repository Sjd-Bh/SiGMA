import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import write
import os
import sys
np.random.seed(42)  # Set a specific seed value, such as 42
import random
random.seed(42)
sys.path.append('/home/bahonar/simulation/SingleCellSim')

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
    end = A_c['startPos'] + (A_c['maxLength'] * A_c['direction'])
    valid_indices = np.where((0 <= A_c['startPos']) & (A_c['startPos'] < refSeq_length) & (0 <= end) & (end < refSeq_length))
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
    A_i['endPos'][mask] = (A_i['startPos'] + A_i['maxLength'] * A_i['direction'])[mask]
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
def updateErrors(A,A_c,parent,end,beta,n_i,refSeq):
    
    # 1. inherit errors from parent and write to amplicons parent (errors which are located in the range of amplicon)
    parent_errors = A[parent]['errors']
    # # if parent_errors is not empty choose the errors that the position of them is between the start and the end of amplicon (between A_c['startPos'] and end)
    # # then add the chosen errors to A_c['errors']
    if len(parent_errors) > 0:
        selected_errors = []
        for start_pos, end_pos in zip(A_c['startPos'], end):
            # Check A_c_direction and updateErr accordingly
            if (A_c['direction'] == -1).any():
                errors_in_range = [(pos, nucleotide) for pos, nucleotide in parent_errors if end_pos <= pos <= start_pos]
            else:
                errors_in_range = [(pos, nucleotide) for pos, nucleotide in parent_errors if start_pos <= pos <= end_pos]        
            # Add the selected errors to the list
            selected_errors.append(errors_in_range)        
        A_c['errors'] = selected_errors
    # # 2. find the numbers of new errors positions in amplicons based on a binomial on the length of amplicon and the error amplification rate beta
    error_positions_count = np.random.binomial(A_c['maxLength'], beta, n_i)
    # # and find the positions of those errors on amplicon randomly
    temp_start = np.minimum(A_c['startPos'], end)
    temp_end = np.maximum(A_c['startPos'], end)   
    selectedPositions = [np.random.randint(start, end + 1, size=count) for start, end, count in zip(temp_start, temp_end, error_positions_count)]
    # # 3. find alternate nucleotide with getAlt for each error position
    A_c['errors'] = [[(pos, getAlt(refSeq[pos])) for pos in positions] for positions in selectedPositions]
    return A_c
####################################################################################
# MDASimulation: this function get the positions of amplicons that is produced in the amplification simulation
# Theta: the number of nucleotides that is made by polymerase in MDA time unit based on polymerase speed
# Gamma: the number of nucleotides that is displaced 
def MDASimulation(patSeq, matSeq, Theta=3000, Gamma=50, DNACoef=3,
                            lMin=2000, lMax=70000, Lambda=0.0001,
                            delta_t=3, beta=0.0001, exclude=False,
                            saveInterval=15, output_folder="output", resume=False,
                            depth = 30):
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
    check_DNA = P
    initial_DNA = check_DNA
    final_DNA = check_DNA * DNACoef
    A = np.array([
        (1, M, M, +1, -1, [], True, 0.0, 'M'),
        (M, 1, M, -1, -1, [], True, 0.0, 'M'),
        (1, P, P, +1, -1, [], True, 0.0, 'P'),
        (P, 1, P, -1, -1, [], True, 0.0, 'P')
    ], dtype=main_dtype)
    At = A
    t = 0
    total_maxLength = 0 
    num_amplicons = 0
    while check_DNA <= final_DNA:
        # Generate new amplicons for all parents using vectorized approach
        n_i_values = np.random.poisson(np.abs(Lambda * (At['endPos'] - At['startPos']) * At['direction']))
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
        num_amplicons += len(newAmplicons)
        # print(t)
        check_DNA += np.sum(newAmplicons['maxLength'])
        
    print(t)
    
    # print(num_amplicons)
    average_maxLength = total_maxLength / num_amplicons
    # print("Average maxLength:", average_maxLength)
    # subseting amplicons for the desired depth 
    subset_percentage = int((depth*initial_DNA)/average_maxLength)
    # print(subset_percentage)
    # num_elements_to_subset = int(len(A) * subset_percentage)
    # print(num_elements_to_subset)
    # print(len(A))
    random_indices = np.random.choice(np.arange(4,len(A)), subset_percentage, replace=False)
    selected_indices = np.concatenate(([0, 1, 2, 3], random_indices))
    subsetA = A[selected_indices]
    # print(len(subsetA))
    # total_length_maxLength = np.sum(subsetA['maxLength'])
    # print(total_length_maxLength)
    return subsetA 
 
####################################################################################
# subsetAmpliconSaveTofFASTA: this function first subset a desired percentage of amplicons randomly
# then convert them to a fasta file that is ready for the next step, sequencing by Art
def subsetAmpliconSaveToFASTA(amplicons, patSeq, matSeq, output_folder="output"):

    records = []
    len_patSeq = len(patSeq)
    len_matSeq = len(matSeq)
    # Add patSeq and its complement
    # records.append(SeqRecord(Seq(patSeq), id="patSeq_1", description=f"Start: 0 , End: {len_patSeq}"))
    # records.append(SeqRecord(Seq(patSeq), id="patSeq_-1", description=f"Start: 0 , End: {len_patSeq}"))

    # # Add matSeq and its complement
    # records.append(SeqRecord(Seq(matSeq), id="matSeq_1", description=f"Start: 0 , End: {len_matSeq}"))
    # records.append(SeqRecord(Seq(matSeq), id="matSeq_-1", description=f"Start: 0 , End: {len_matSeq}"))

    for amplicon in amplicons:
        length = amplicon['maxLength']
        direction = amplicon['direction']

        if direction == -1:
            startPos  = amplicon['startPos'] + length*direction + 1 
            endPos = amplicon['startPos'] + 1
        else:
            startPos = amplicon['startPos']
            endPos = amplicon['startPos'] + length*direction + 1
                    
        refSeq = patSeq if amplicon['source'] == 'P' else matSeq
        seqName = "patSeq" if amplicon['source'] == 'P' else "matSeq"
        
        # ampliconSeq = refSeq[startPos:endPos + direction:direction]
        # ampliconSeq = applyErrors(ampliconSeq, amplicon['errors'], direction, length, startPos,endPos )
        ampliconSeq = refSeq[startPos:endPos]

        record_id = f"{seqName}_{direction}"
        record = SeqRecord(Seq(ampliconSeq), id=record_id, description=f"Start: {startPos}, End: {endPos - 1}")

        records.append(record)

    outputFilename = os.path.join(output_folder, "subset.fasta")
    write(records, outputFilename, "fasta")

#####################################################################################