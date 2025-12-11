import pandas as pd
import pysam
from collections import Counter

def analyze_editing_events(bam_file, chrom, start, end, meta_file=None, status=None, window=3):
    """
    Analyze deduplicated reads to detect editing events (KO/KI) using CIGAR.
    Extended to:
      1) count number of edited bases per read
      2) compute proportion of edit types across edited reads
      3) compute detailed counts of base-length-specific events (I1, I2, D1, S3...)
    """

    start_ext = max(0, start - window)
    end_ext = end + window

    barcodes_to_use = None
    if meta_file:
        meta_df = pd.read_csv(meta_file)
        barcodes_to_use = set(meta_df.loc[meta_df['status'] == status, 'barcode'])

    bam = pysam.AlignmentFile(bam_file, "rb")
    umi_dict = {}

    # For global editing proportions
    edit_type_counter = Counter()          # counts of I, D, S
    edit_length_counter = Counter()        # counts of I_1, I_2, D_3, S_5...

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
        if key not in umi_dict:
            cigar_info = []
            ref_pos = read.reference_start

            edited_bases = 0      # total edited bases (I + D + S)
            edit_types = []       # record types

            for op, length in read.cigartuples:

                if op == 0:  # M
                    cigar_info.append(('M', ref_pos, ref_pos + length))
                    ref_pos += length

                elif op == 1:  # I
                    cigar_info.append(('I', ref_pos, length))
                    edited_bases += length
                    edit_types.append('I')

                    edit_length_counter[f"I_{length}"] += 1
                    edit_type_counter["I"] += 1

                elif op == 2:  # D
                    cigar_info.append(('D', ref_pos, ref_pos + length))
                    edited_bases += length
                    edit_types.append('D')

                    edit_length_counter[f"D_{length}"] += 1
                    edit_type_counter["D"] += 1

                    ref_pos += length

                elif op == 4:  # S
                    cigar_info.append(('S', ref_pos, length))
                    edited_bases += length
                    edit_types.append('S')

                    edit_length_counter[f"S_{length}"] += 1
                    edit_type_counter["S"] += 1

                elif op == 3:  # N
                    cigar_info.append(('N', ref_pos, ref_pos + length))
                    ref_pos += length

                else:
                    ref_pos += length

            umi_dict[key] = {
                'name': read.query_name,
                'start': read.reference_start,
                'end': read.reference_end,
                'cigar': cigar_info,
                'edited_bases': edited_bases,
                'edit_types': edit_types
            }

    bam.close()

    # ---- compute proportions ----
    total_edits = sum(edit_type_counter.values())
    if total_edits > 0:
        edit_proportions = {
            etype: count / total_edits
            for etype, count in edit_type_counter.items()
        }
    else:
        edit_proportions = {"I": 0, "D": 0, "S": 0}

    return list(umi_dict.values()), edit_proportions, edit_length_counter
