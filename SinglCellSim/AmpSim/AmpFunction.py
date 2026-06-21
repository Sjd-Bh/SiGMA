import os
import pickle
import csv
import random
import numpy as np
from Bio.Seq import Seq
import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter
from Bio import SeqIO

####################################################
# CNV Loading and Breakpoint Arrays
####################################################
def load_cnvs_from_bed(bed_file):
    """
    Load CNVs from BED-like file and assign them to paternal/maternal alleles.
    Returns:
        cnvs: full CNV info for each allele
        cnvs_for_mda: minimal info for MDA simulation
    """
    cnvs = []
    with open(bed_file, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split() 
            chrom, start, end, cnv_id, cnv_type, allele = parts[:6]

    	    # Hardcode prevalence to 1.0 (100%) since it's not in the BED file
            prevalence = 1.0
            start, end = int(start), int(end)
            #prevalence = float(prev)
            # Assign alleles based on prevalence
            alleles = ["P", "M"] if prevalence == 1 else ["P"] if random.random() < prevalence else ["M"]
            for allele in alleles:
                if cnv_type == "DEL":
                    cnvs.append({"type": "DEL", "start": start, "end": end, "allele": allele, "copy_number": 0})
                elif cnv_type == "DUP":
                    copies = max(2, int(round(prevalence * 2)))
                    cnvs.append({"type": "DUP", "start": start, "end": end, "allele": allele, "copy_number": copies})
    cnvs_for_mda = [{"start": c["start"], "end": c["end"], "copy_number": c.get("copy_number", 1), "allele": c.get("allele","P")} for c in cnvs]
    return cnvs, cnvs_for_mda

def build_breakpoints(cnv_regions, genome_length):
    """
    Converts CNV regions into a Breakpoint Array.
    Returns:
        breakpoints: np.array of positions where copy number changes
        cn_states: np.array of copy numbers corresponding to the segments between breakpoints
        dense_cn: np.array of length genome_length containing per-base copy number
    """
    dense_cn = np.ones(genome_length, dtype=int)
    for c in cnv_regions:
        s, e, cn = c['start'], c['end'], c.get('copy_number', 1)
        # Apply CNV (bounds checking)
        s, e = max(0, s), min(genome_length, e)
        dense_cn[s:e] = cn
        
    # Find points where copy number changes
    changes = np.where(dense_cn[:-1] != dense_cn[1:])[0] + 1
    
    # Construct breakpoint arrays
    breakpoints = np.concatenate(([0], changes, [genome_length]))
    cn_states = np.concatenate(([dense_cn[0]], dense_cn[changes], [0]))
    
    return breakpoints, cn_states, dense_cn

def get_copy_number(positions, breakpoints, cn_states):
    """ Fast O(log N) lookup of copy number at reference positions """
    indices = np.searchsorted(breakpoints, positions, side='right') - 1
    return cn_states[indices]

def get_extension_limit(pos, direction, breakpoints, cn_states):
    """ 
    Finds the maximum reference coordinate an amplicon can reach 
    before hitting a homozygous deletion boundary or genome edge.
    """
    idx = np.searchsorted(breakpoints, pos, side='right') - 1
    if cn_states[idx] == 0:
        return pos  # Already in a deletion, cannot extend
        
    deletion_indices = np.where(cn_states == 0)[0]
    
    if direction == 1:
        # Find the start of the next deletion
        del_starts = breakpoints[deletion_indices]
        valid = del_starts > pos
        return del_starts[valid][0] if np.any(valid) else breakpoints[-1]
    else:
        # Find the end of the previous deletion
        # FIX: Ignore the very last index since it's the genome boundary and has no right-side breakpoint
        safe_indices = deletion_indices[deletion_indices < len(breakpoints) - 1]
        del_ends = breakpoints[safe_indices + 1] 
        valid = del_ends <= pos
        return del_ends[valid][-1] if np.any(valid) else 0


####################################################
# Amplicon Utilities
####################################################
def amp_span(amp):
    s = min(amp['startPos'], amp['endPos'])
    e = max(amp['startPos'], amp['endPos'])
    return s, e, e - s

def getAlt(base):
    return {'A':'C','C':'T','G':'A','T':'G','N':'N'}.get(base, 'N')

# --------------------------------------------------
# Update Polya urn counts
# --------------------------------------------------
def update_polya_counts(A, P_count, M_count, lMin):
    """
    Update dense arrays natively in reference coordinates.
    """
    for amp in A:
        if amp['released']:
            continue
        s, e, length = amp_span(amp)
        if length < lMin:
            continue

        if amp['source'] == 'P':
            P_count[s:e] += 1
        else:
            M_count[s:e] += 1
            
    return M_count, P_count

####################################################
# Amplicon Extension
####################################################
def extendAmplicon(amplicons, delta_t, Theta, bp_P, cn_P, bp_M, cn_M, beta, refSeq_P, refSeq_M, return_extended=False):
    """ Extend amplicons directly in reference coordinates until limits or boundaries, introducing errors """
    increment = int(np.floor(Theta * delta_t))
    if increment <= 0 or len(amplicons) == 0:
        return (amplicons, np.array([], dtype=amplicons.dtype)) if return_extended else amplicons

    new_regions = []

    for amp in amplicons:
        if amp['released']:
            continue

        bp = bp_P if amp['source'] == 'P' else bp_M
        cn = cn_P if amp['source'] == 'P' else cn_M
        refSeq = refSeq_P if amp['source'] == 'P' else refSeq_M

        cur_len = abs(amp['endPos'] - amp['startPos'])
        new_len = min(cur_len + increment, amp['maxLength'])
        delta = new_len - cur_len
        if delta <= 0:
            continue

        old_end = amp['endPos']
        limit = get_extension_limit(old_end, amp['direction'], bp, cn)
        
        # Extend considering deletion boundaries
        if amp['direction'] == 1:
            amp['endPos'] = min(old_end + delta, limit)
            lo, hi = old_end, amp['endPos'] - 1
        else:
            amp['endPos'] = max(old_end - delta, limit)
            lo, hi = amp['endPos'], old_end - 1

        # Calculate actual extension length to generate errors
        actual_delta = abs(amp['endPos'] - old_end)
        
        # --- Introduce New Errors ---
        if actual_delta > 0 and beta > 0:
            num_errors = np.random.binomial(actual_delta, beta)
            if num_errors > 0:
                error_positions = np.random.randint(lo, hi + 1, num_errors)
                new_errors = []
                for pos in error_positions:
                    alt_base = getAlt(refSeq[pos])
                    new_errors.append((pos, alt_base))
                
                # Append to the amplicon's existing error list
                amp['errors'].extend(new_errors)

        # Mark released if boundary or maxLength reached
        if abs(amp['endPos'] - amp['startPos']) >= amp['maxLength'] or amp['endPos'] == limit:
            amp['released'] = True

        if return_extended:
            ext_start, ext_end = (old_end, amp['endPos']) if amp['direction'] == 1 else (amp['endPos'], old_end)
            new_regions.append((amp['startPos'], amp['endPos'],
                                abs(ext_end - ext_start), amp['direction'], amp['parent'],
                                amp['errors'], amp['released'], amp['startTime'],
                                amp['source'], amp['copy_index']))

    if return_extended:
        A_ext = np.array(new_regions, dtype=amplicons.dtype) if new_regions else np.array([], dtype=amplicons.dtype)
        return amplicons, A_ext

    return amplicons


####################################################
# Error Update
####################################################
def updateErrors(A, A_c, parent_indices):
    """ Repurposed to STRICTLY handle the inheritance of parental errors. """
    for i in range(len(A_c)):
        p_idx = parent_indices[i]
        p_errors = A[p_idx]['errors']
        
        if p_errors:
            start = A_c['startPos'][i]
            direction = A_c['direction'][i]
            max_l = A_c['maxLength'][i]
            
            # Filter parent errors to only those within the child's potential extension range
            if direction == 1:
                lo, hi = start, start + max_l
            else:
                lo, hi = start - max_l, start
                
            inherited = [(pos, alt) for (pos, alt) in p_errors if lo <= pos <= hi]
            A_c['errors'][i] = inherited
        else:
            A_c['errors'][i] = []
            
    return A_c

####################################################
# Generate New Amplicons
####################################################
def GenerateNewAmp(
    refSeq_P, refSeq_M, A, parent_templates, t, delta_t, main_dtype,
    lengths, positions, parent_indices, beta, Theta, bp_P, cn_P,
    bp_M, cn_M, genome_length, copy_indices
):
    n_events = len(positions)
    if n_events == 0:
        return np.array([], dtype=main_dtype)

    A_c = np.zeros(n_events, dtype=main_dtype)
    A_c['source']     = parent_templates['source']
    A_c['direction']  = -parent_templates['direction']
    A_c['parent']     = parent_indices

    parent_start = A[parent_indices]['startPos']
    distance_from_parent = np.abs(parent_start - positions)
    
    A_c['startPos'] = positions
    A_c['endPos'] = positions
    
    limited_lengths = np.minimum(lengths, distance_from_parent)
    
    for i in range(n_events):
        bp = bp_P if A_c['source'][i] == 'P' else bp_M
        cn = cn_P if A_c['source'][i] == 'P' else cn_M
        limit = get_extension_limit(positions[i], A_c['direction'][i], bp, cn)
        dist_to_limit = abs(limit - positions[i])
        limited_lengths[i] = min(limited_lengths[i], dist_to_limit)

    A_c['maxLength'] = limited_lengths

    valid_mask = limited_lengths >= 250
    if not np.any(valid_mask):
        return np.array([], dtype=main_dtype)

    A_c = A_c[valid_mask]
    parent_indices  = parent_indices[valid_mask]
    copy_indices    = copy_indices[valid_mask]
    n_events = len(A_c)
    
    A_c['startTime']  = np.random.uniform(t, t + delta_t, n_events)
    A_c['released']   = False
    A_c['copy_index'] = copy_indices
    
    empty_errors = np.empty(n_events, dtype=object)
    for i in range(n_events): empty_errors[i] = []
    A_c['errors'] = empty_errors

    # Inherit errors directly from parents (No new errors generated here)
    valid = ~A_c['released']
    if np.any(valid):
        A_c[valid] = updateErrors(A, A_c[valid], parent_indices[valid])

    return A_c

def apply_cnvs_to_sequence(seq, cnvs, seq_id="genome"):
    """ Simple placeholder for sequence handling if needed outside simulation """
    if isinstance(seq, SeqRecord): ref_seq = str(seq.seq)
    else: ref_seq = str(seq)
    # With purely reference simulation, we don't strictly generate a CNV modified string here anymore,
    # as FASTA output applies the true reference and errors natively.
    return SeqRecord(Seq(ref_seq), id=seq_id, description="Reference")

####################################################
# FASTA Output
####################################################
def subsetAmpliconSaveToFASTA(amplicons, refSeq_P, refSeq_M, lMin=500, output_folder="output"):
    """
    Directly slice from reference using startPos and endPos and apply track errors.
    No effective coordinate conversion needed.
    """
    os.makedirs(output_folder, exist_ok=True)
    out_fasta = os.path.join(output_folder, "subset.fasta")
    records = []
    
    # Filter roots and short sequences
    valid_amps = amplicons[(amplicons['parent'] != -1) & 
                           (np.abs(amplicons['endPos'] - amplicons['startPos']) >= lMin)]
    
    for i, amp in enumerate(valid_amps):
        start = min(amp['startPos'], amp['endPos'])
        end = max(amp['startPos'], amp['endPos'])
        
        seq_list = list(refSeq_P[start:end] if amp['source'] == 'P' else refSeq_M[start:end])
   
        if amp['errors']:
            for err_pos, alt_base in amp['errors']:
                local_idx = err_pos - start
                if 0 <= local_idx < len(seq_list):
                    seq_list[local_idx] = alt_base
                    
        final_seq = "".join(seq_list)
        
        if amp['direction'] == -1:
            final_seq = str(Seq(final_seq).reverse_complement())
            
        allele = "pat" if amp['source'] == "P" else "mat"
        record_id = f"Amp_{i}_Src:{allele}_Copy:{amp['copy_index']}_Pos:{start}-{end}_Dir:{amp['direction']}"
        records.append(SeqRecord(Seq(final_seq), id=record_id, description=""))
        
    with open(out_fasta, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")
        
    print(f"Saved {len(records)} amplicons to {out_fasta}")
    return out_fasta

####################################################
# Parent Template Selection
####################################################
def find_template_index(A, ref_positions, sources, copy_indices, min_dist=2500, max_dist=70000, min_length=1500):
    n = len(ref_positions)
    parent_indices = np.full(n, -1, dtype=int)

    if len(A) == 0 or n == 0:
        return parent_indices

    A_start = np.minimum(A['startPos'], A['endPos'])
    A_end   = np.maximum(A['startPos'], A['endPos'])
    A_len = np.abs(A['endPos'] - A['startPos'])

    for i, pos in enumerate(ref_positions):
        allele = sources[i]
        copy   = copy_indices[i]

        mask = (
            (A['source'] == allele) &
            ((A['copy_index'] == copy) | (A['copy_index'] == -1)) &
            (A_start <= pos) &
            (A_end >= pos) &
            (A_len > min_length)
        )

        if not np.any(mask): continue
        candidates = np.where(mask)[0]

        dist = np.abs(A['startPos'][candidates] - pos)
        ok = dist >= min_dist

        candidates = candidates[ok]
        if len(candidates) > 0:
            parent_indices[i] = np.random.choice(candidates)
    
    return parent_indices

def collapse_copycount_regions(P_count, M_count):
    P_count = np.asarray(P_count)
    M_count = np.asarray(M_count)
    genome_length = len(P_count)

    change = np.zeros(genome_length, dtype=bool)
    change[0] = True
    change[1:] = (P_count[1:] != P_count[:-1]) | (M_count[1:] != M_count[:-1])

    starts = np.where(change)[0]
    ends   = np.r_[starts[1:], genome_length]

    rows = []
    for s, e in zip(starts, ends):
        rows.append({"start": s, "end": e - 1, "P_count": int(P_count[s]), "M_count": int(M_count[s]), "length": e - s})
    return pd.DataFrame(rows)


###################################################
# MDA Simulation
####################################################
def MDASimulation(
    patSeq, matSeq, lMin, lMax, delta_t, beta, CNVs=None,
    output_folder="output", total_time=10.0, t_switch=2.0,
    lambda_init=1e-7, lambda_exp=1e-6, Theta_init=50, Theta_exp=500,
    main_dtype=None
):
    genome_length = len(patSeq)
    CNVs = [] if CNVs is None else CNVs
    os.makedirs(output_folder, exist_ok=True)

    if main_dtype is None:
        main_dtype = np.dtype([
            ('startPos', int), ('endPos', int),
            ('maxLength', int), ('direction', int),
            ('parent', int), ('errors', 'O'),
            ('released', bool), ('startTime', float),
            ('source', 'U1'), ('copy_index', int)
        ])

    cnvs_P = [c for c in CNVs if c.get("allele","P")=="P"]
    cnvs_M = [c for c in CNVs if c.get("allele","M")=="M"]
    
    bp_P, cn_P, P_count = build_breakpoints(cnvs_P, genome_length)
    bp_M, cn_M, M_count = build_breakpoints(cnvs_M, genome_length)

    A = np.array([
        (0, genome_length-1, genome_length, +1, -1, [], True, 0.0, 'P', -1),
        (genome_length-1, 0, genome_length, -1, -1, [], True, 0.0, 'P', -1),
        (0, genome_length-1, genome_length, +1, -1, [], True, 0.0, 'M', -1),
        (genome_length-1, 0, genome_length, -1, -1, [], True, 0.0, 'M', -1)
    ], dtype=main_dtype)

    t = 0.0
    cycle_idx = 0
    cycle_stats = []
    archived_amplicons_per_cycle = []

    while t < total_time:
        if t < t_switch:
            phase = "initiation"
            lambda_rate = lambda_init
            Theta = Theta_init
            template_pool = A
        else:
            phase = "exponential"
            lambda_rate = lambda_exp
            Theta = Theta_exp
            non_root_amplicons = A[A["parent"] != -1]
            template_pool = non_root_amplicons[non_root_amplicons["maxLength"] > 1000]
            
        if len(template_pool) == 0:
            archived_amplicons_per_cycle.append(np.array([], dtype=main_dtype))
            cycle_stats.append((cycle_idx, phase, 0, 0))
            t += delta_t
            cycle_idx += 1
            continue

        n_i = np.random.poisson(lambda_rate * genome_length * delta_t)
        if n_i == 0:
            archived_amplicons_per_cycle.append(np.array([], dtype=main_dtype))
            cycle_stats.append((cycle_idx, phase, 0, 0))
            t += delta_t
            cycle_idx += 1
            continue

        if phase == "initiation":
            ref_positions = np.random.randint(0, genome_length, n_i)
        else:
            chosen_idx = np.random.randint(0, len(template_pool), n_i)
            ref_positions = np.zeros(n_i, dtype=int)
            for i, idx in enumerate(chosen_idx):
                start, end = template_pool[idx]['startPos'], template_pool[idx]['endPos']
                
                # Prevent upper bound from exceeding genome_length - 1
                upper_bound = min(max(start, end), genome_length - 1) 
                ref_positions[i] = np.random.randint(min(start, end), upper_bound + 1) 
        
        ref_positions = np.clip(ref_positions, 0, genome_length - 1)
        denom = P_count[ref_positions] + M_count[ref_positions]
        p_prob = np.divide(P_count[ref_positions], denom, out=np.zeros_like(denom, dtype=float), where=denom > 0)
        
        rand_vals = np.random.rand(n_i)
        sources = np.where(rand_vals < p_prob, 'P', 'M').astype('<U1')
        
        valid_mask = denom > 0
        sources = sources[valid_mask]
        ref_positions = ref_positions[valid_mask]
        n_i = len(ref_positions)

        copy_indices = np.zeros(n_i, dtype=int)
        for i in range(n_i):
            local_cn = P_count[ref_positions[i]] if sources[i] == 'P' else M_count[ref_positions[i]]
            copy_indices[i] = np.random.randint(1, local_cn + 1) if local_cn > 0 else -1

        if phase == "initiation" and lMax < 2000:
            lengths = np.random.randint(2000, 10000 + 1, n_i)
        else:
            lengths = np.random.randint(lMin, lMax + 1, n_i)

        parent_indices = find_template_index(template_pool, ref_positions, sources, copy_indices)
        valid = parent_indices != -1
        A_c = np.array([], dtype=main_dtype)

        if np.any(valid):
            A_c = GenerateNewAmp(
                refSeq_P=patSeq, refSeq_M=matSeq, A=A,
                parent_templates=template_pool[parent_indices[valid]],
                t=t, delta_t=delta_t, main_dtype=main_dtype,
                lengths=lengths[valid], positions=ref_positions[valid],
                parent_indices=parent_indices[valid], beta=beta, Theta=Theta,
                bp_P=bp_P, cn_P=cn_P, bp_M=bp_M, cn_M=cn_M,
                genome_length=genome_length, copy_indices=copy_indices[valid]
            )
            
            if len(A_c) > 0:
                A = np.concatenate((A, A_c))
                archived_amplicons_per_cycle.append(A_c.copy())
            else:
                archived_amplicons_per_cycle.append(np.array([], dtype=main_dtype))

        # UPDATED extendAmplicon call with beta, patSeq, and matSeq
        A, A_ext = extendAmplicon(
            A, delta_t, Theta, bp_P, cn_P, bp_M, cn_M, 
            beta, patSeq, matSeq, return_extended=True
        )

        if len(A_c) > 0 or len(A_ext) > 0:
            combined = np.concatenate([A_c, A_ext]) if len(A_c) > 0 else A_ext
            M_count, P_count = update_polya_counts(combined, P_count, M_count, lMin)
            
        total_len = (sum(A_c['maxLength']) if len(A_c) > 0 else 0) + (sum(A_ext['maxLength']) if len(A_ext) > 0 else 0)
        cycle_stats.append((cycle_idx, phase, len(A_c), total_len))

        t += delta_t
        cycle_idx += 1

    # Add this at the very end of MDASimulation to count raw mathematical errors
    total_absolute_errors = sum(len(amp['errors']) for amp in A)
    print(f"Total raw error events simulated: {total_absolute_errors}")


    with open(os.path.join(output_folder, "amplicon_stats.tsv"), "w", newline='') as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["Cycle", "Phase", "NewAmplicons", "TotalLength_bp"])
        writer.writerows(cycle_stats)

    with open(os.path.join(output_folder, "amplicons.pkl"), "wb") as f:
        pickle.dump(A, f)

    return A, cycle_stats, P_count, M_count, CNVs, bp_P, bp_M
