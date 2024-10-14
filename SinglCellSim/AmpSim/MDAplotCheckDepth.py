import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import write
import os
import matplotlib.pyplot as plt
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
def MDASimulation(patSeq, matSeq, Theta=12000, Gamma=50, DNACoef=400,
                  lMin=2000, lMax=70000, Lambda=0.0001,
                  delta_t=0.01, beta=0.000001, exclude=False,
                  saveInterval=15, output_folder="output", resume=False,
                  amp_depth = 5, template = False):
    P = len(patSeq)
    M = len(matSeq)
    genome_length = P  # Assuming patSeq and matSeq are of the same length
    coverage = np.zeros(genome_length)  # Initialize coverage array to track depth
    coverage_cycle = []  # List to store coverage at each cycle
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    # Other parts of the function remain the same...   
    main_dtype = np.dtype([('startPos', int),
                           ('endPos', int),
                           ('maxLength', int),
                           ('direction', int),
                           ('parent', int),
                           ('errors', 'O'),
                           ('released', bool),
                           ('startTime', float),
                           ('source', 'U1')])
    check_DNA = 0
    initial_DNA = P
    final_DNA = P * DNACoef
    A = np.array([(0, M-1, M, +1, -1, [], True, 0.0, 'M'),
                  (M-1, 1, M, -1, -1, [], True, 0.0, 'M'),
                  (0, P-1, P, +1, -1, [], True, 0.0, 'P'),
                  (P-1, 0, P, -1, -1, [], True, 0.0, 'P')],
                 dtype=main_dtype)
    t = 0
    total_maxLength = 0 
    num_amplicons = 0
    At = A
    cycle_count = 0    
    while check_DNA <= final_DNA:
        n_i_values = np.random.poisson(np.abs(Lambda * delta_t * (At['endPos'] - At['startPos']) * At['direction']))
        newAmplicons = np.concatenate([GenerateNewAmp(patSeq if At[parent]['source'] == 'P' else matSeq,
                                                      At, At[parent], t, delta_t, main_dtype, lMin, lMax, parent, Gamma, beta, Theta, exclude, n_i)
                                       for parent, n_i in enumerate(n_i_values)])
        unreleased_mask = ~A['released']
        A[unreleased_mask] = extendAmplicon(A[unreleased_mask], t, delta_t, Theta)        
        # Update coverage array based on amplicon start and end positions
        for amp in newAmplicons:
            coverage[amp['startPos']:amp['endPos']] += 1        
        A = np.concatenate((A, newAmplicons))
        t += delta_t
        total_maxLength += np.sum(newAmplicons['maxLength'])
        num_amplicons += len(newAmplicons)
        check_DNA += np.sum(newAmplicons['maxLength'])        
        # Save coverage and plot after each cycle
        if cycle_count % saveInterval == 0:
            coverage_cycle.append(coverage.copy())  # Store a snapshot of the coverage
            plt.figure()
            plt.plot(coverage)
            plt.xlabel('Genome Positions')
            plt.ylabel('Depth of Coverage')
            plt.title(f'Depth of Coverage - Cycle {cycle_count}')
            plt.savefig(os.path.join(output_folder, f"coverage_cycle_{cycle_count}.png"))
            plt.close()  # Close the figure to avoid displaying it during the run
        cycle_count += 1    
    # Save the final state and plot
    plt.plot(coverage)
    plt.xlabel('Genome Positions')
    plt.ylabel('Depth of Coverage')
    plt.title('Final Depth of Coverage Across Genome Positions')
    final_plot_path = os.path.join(output_folder, "final_coverage.png")
    plt.savefig(final_plot_path)
    plt.close()  # Close the final figure
    # Save the amplicons to a pickle file
    # outputFilename = os.path.join(output_folder, "amplicons.pkl")
    # with open(outputFilename, 'wb') as f:  
    #     pickle.dump(A, f)
    return A

####################################################################################
patSeq = "../CellGenome/test2_afterEditingZero_basedSys/200/paternal_cell.fasta"
matSeq = "../CellGenome/test2_afterEditingZero_basedSys/200/maternal_cell.fasta"
MDASimulation(patSeq, matSeq, Theta=12000, Gamma=50, DNACoef=5,lMin=2000, lMax=70000,
              Lambda=0.0001,delta_t=0.01, beta=0.000001, exclude=False,
              saveInterval=15, output_folder="../test/plotsDP2", resume=False,
              amp_depth=5, template=True)
