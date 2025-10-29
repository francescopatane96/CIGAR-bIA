import pandas as pd
import pysam

def analyze_editing_events(bam_file, chrom, start, end, meta_file=None, status=None, window=3):
    """
    Analyze deduplicated reads to detect editing events (KO/KI) using CIGAR.

    Args:
        bam_file (str): Path to indexed BAM file.
        chrom (str): Chromosome.
        start (int): Start position of the target region.
        end (int): End position of the target region.
        meta_file (str, optional): CSV file with 'barcode' and 'status' columns.
        status (str, optional): Status to filter in meta_file.
        window (int): Extension of the window around start-end.

    Returns:
        parsed_reads (list of dict): List of reads with CIGAR info, deduplicated by UMI.
    """
    start_ext = max(0, start - window)
    end_ext = end + window

    barcodes_to_use = None
    if meta_file:
        meta_df = pd.read_csv(meta_file)
        barcodes_to_use = set(meta_df.loc[meta_df['status'] == status, 'barcode'])

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
        if key not in umi_dict:
            cigar_info = []
            ref_pos = read.reference_start
            for op, length in read.cigartuples:
                if op == 0:  # M
                    cigar_info.append(('M', ref_pos, ref_pos + length))
                    ref_pos += length
                elif op == 1:  # I
                    cigar_info.append(('I', ref_pos, length))
                elif op == 2:  # D
                    cigar_info.append(('D', ref_pos, ref_pos + length))
                    ref_pos += length
                elif op == 3:  # N
                    cigar_info.append(('N', ref_pos, ref_pos + length))
                    ref_pos += length
                elif op == 4:  # S
                    cigar_info.append(('S', ref_pos, length))
                else:
                    ref_pos += length

            umi_dict[key] = {
                'name': read.query_name,
                'start': read.reference_start,
                'end': read.reference_end,
                'cigar': cigar_info
            }

    bam.close()
    return list(umi_dict.values())
