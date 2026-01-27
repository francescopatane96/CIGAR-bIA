import pandas as pd
import pysam
from collections import Counter

def analyze_editing_events(bam_file, chrom, start, end, meta_file=None, status_col=None, status=None, window=3):
    """
    Analizza reads deduplicate per (barcode, UMI) e rileva editing events.
    Ora considera I, D e N come modifiche.
    
    Restituisce lista di dict con:
        - name, start, end
        - cigar
        - I_count, D_count, N_count
        - edited (True/False)
        - edited_bases
        - frameshift (True se non multiplo di 3)
        - barcode, umi
    """
    start_ext = max(0, start - window)
    end_ext = end + window

    barcodes_to_use = None
    if meta_file:
        meta_df = pd.read_csv(meta_file)
        barcodes_to_use = set(meta_df.loc[meta_df[status_col] == status, 'barcode'])

    bam = pysam.AlignmentFile(bam_file, "rb")
    umi_dict = {}

    for read in bam.fetch(chrom, start_ext, end_ext):
        if read.is_unmapped or read.cigartuples is None:
            continue
        try:
            bc = read.get_tag('CB')
            umi = read.get_tag('UB')
        except KeyError:
            continue

        if barcodes_to_use and bc not in barcodes_to_use:
            continue

        key = (bc, umi)
        if key in umi_dict:
            continue  # deduplicazione

        cigar_info = []
        ref_pos = read.reference_start
        edited_bases = 0
        I_count = D_count = N_count = 0

        for op, length in read.cigartuples:
            if op == 0:  # M
                cigar_info.append(('M', ref_pos, ref_pos + length))
                ref_pos += length
            elif op == 1:  # I
                cigar_info.append(('I', ref_pos, length))
                edited_bases += length
                I_count += length
            elif op == 2:  # D
                cigar_info.append(('D', ref_pos, ref_pos + length))
                edited_bases += length
                D_count += length
                ref_pos += length
            elif op == 3:  # N (skipped)
                cigar_info.append(('D', ref_pos, ref_pos + length))
                edited_bases += length
                N_count += length
                ref_pos += length
            elif op == 4:  # S
                cigar_info.append(('S', ref_pos, length))
                ref_pos += length
            else:
                ref_pos += length

        frameshift = (edited_bases % 3) != 0
        edited = edited_bases > 0

        umi_dict[key] = {
            'name': read.query_name,
            'start': read.reference_start,
            'end': read.reference_end,
            'cigar': cigar_info,
            'I_count': I_count,
            'D_count': D_count,
            'N_count': N_count,
            'edited_bases': edited_bases,
            'edited': edited,
            'frameshift': frameshift,
            'barcode': bc,
            'umi': umi
        }

    bam.close()
    return list(umi_dict.values())
