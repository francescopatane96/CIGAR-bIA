import pandas as pd
import pysam
from collections import Counter

def analyze_and_save_editing_events(
    bam_file,
    chrom,
    start,
    end,
    output_csv,
    meta_file=None,
    status=None,
    window=3
):
    """
    Analizza tutte le reads in una regione, salva in un unico CSV:
      - se la read è editata
      - numero totale di basi editate
      - conteggio I/D/S
      - se provoca frameshift
    """

    start_ext = max(0, start - window)
    end_ext = end + window

    barcodes_to_use = None
    if meta_file:
        meta = pd.read_csv(meta_file)
        barcodes_to_use = set(meta.loc[meta["status"] == status, "barcode"])

    bam = pysam.AlignmentFile(bam_file, "rb")

    results = []

    for read in bam.fetch(chrom, start_ext, end_ext):

        if read.is_unmapped or read.cigartuples is None:
            continue

        try:
            bc = read.get_tag("CB")
            umi = read.get_tag("UB")
        except KeyError:
            continue

        if barcodes_to_use and bc not in barcodes_to_use:
            continue

        # ---- ANALISI CIGAR ----
        edit_counter = Counter()
        edited_bases = 0

        for op, length in read.cigartuples:

            if op == 1:       # Insertion
                edit_counter["I"] += 1
                edited_bases += length

            elif op == 2:     # Deletion
                edit_counter["D"] += 1
                edited_bases += length

            elif op == 4:     # Soft clipping
                edit_counter["S"] += 1
                edited_bases += length

        is_edited = edited_bases > 0
        frameshift = False if edited_bases == 0 else (edited_bases % 3 != 0)

        results.append({
            "barcode": bc,
            "umi": umi,
            "read_name": read.query_name,
            "edited": is_edited,
            "edited_bases": edited_bases,
            "I_count": edit_counter["I"],
            "D_count": edit_counter["D"],
            "S_count": edit_counter["S"],
            "frameshift": frameshift
        })

    bam.close()

    # ---- Salva tutto in CSV ----
    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False)

    return df
