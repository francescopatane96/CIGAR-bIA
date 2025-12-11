
import pandas as pd
import pysam
from collections import Counter

def analyze_editing_events(
    bam_file,
    chrom,
    start,
    end,
    meta_file=None,
    status=None,
    window=3
):
    """
    Analizza editing basato sul CIGAR.
    
    Estrae:
      - numero totale di reads editate
      - tipo di modifica (conteggio I/D/S per read)
      - numero di basi modificate per read (somma)
      - numero di reads con lunghezza modificata multipla di 3 o no
    """

    start_ext = max(0, start - window)
    end_ext = end + window

    # filtra barcodes se necessario
    barcodes_to_use = None
    if meta_file is not None:
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

        # ignora reads non editate
        if edited_bases == 0:
            continue

        # registra la read
        results.append({
            "barcode": bc,
            "umi": umi,
            "read_name": read.query_name,
            "I_count": edit_counter["I"],
            "D_count": edit_counter["D"],
            "S_count": edit_counter["S"],
            "edited_bases": edited_bases,
            "multiple_of_3": (edited_bases % 3 == 0)
        })

    bam.close()

    df = pd.DataFrame(results)

    # ---- METRICHE FINALI ----
    total_edited_reads = len(df)
    reads_mult_3 = df["multiple_of_3"].sum()
    reads_not_mult_3 = total_edited_reads - reads_mult_3

    summary = {
        "total_edited_reads": total_edited_reads,
        "reads_multiple_of_3": reads_mult_3,
        "reads_not_multiple_of_3": reads_not_mult_3
    }

    return df, summary
