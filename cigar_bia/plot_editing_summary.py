import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_editing_summary(parsed_reads, chrom, start, end, status=None):
    """
    Deduplica le reads e mostra grafici di editing:
    1) Efficienza KO (reads editate / tot)
    2) Stacked barplot: distribuzione tipi di editing (I, D, S) e basi modificate
       con colore che indica frameshift o no.
    
    Args:
        parsed_reads (list of dict or DataFrame): output di analyze_and_save_editing_events
        chrom (str): cromosoma
        start (int): inizio regione
        end (int): fine regione
        status (str, optional): nome condizione
    """

    # Converti DataFrame in lista di dict
    if isinstance(parsed_reads, pd.DataFrame):
        parsed_reads = parsed_reads.to_dict(orient='records')

    # Deduplicazione
    seen = set()
    dedup_reads = []
    for read in parsed_reads:
        key = (read.get('barcode'), read.get('umi'))
        if key not in seen:
            seen.add(key)
            dedup_reads.append(read)

    # Trasformiamo in DataFrame per comodità
    df = pd.DataFrame(dedup_reads)

    # --- 1. Efficienza KO ---
    total_reads = len(df)
    edited_reads = df['edited'].sum()
    efficiency = edited_reads / total_reads if total_reads > 0 else 0

    print(f"Chromosome {chrom}:{start}-{end}")
    print(f"Total reads: {total_reads}")
    print(f"Edited reads: {edited_reads}")
    print(f"KO efficiency: {efficiency:.3f}")

    # --- 2. Prepara dati per stacked barplot ---
    # Trasforma in formato lungo
    df_long = df.melt(
        id_vars=['read_name', 'frameshift', 'edited'],
        value_vars=['I_count', 'D_count', 'S_count'],
        var_name='edit_type',
        value_name='count'
    )

    # Filtra solo reads editate
    df_long = df_long[df_long['edited'] == True]

    # Se non ci sono reads editate, fermati
    if df_long.empty:
        print("No edited reads in the specified region.")
        return

    # Stacked barplot: basi modificate per tipo
    plt.figure(figsize=(8, 5))
    sns.barplot(
        data=df_long,
        x='read_name',
        y='count',
        hue='edit_type',
        dodge=False
    )
    plt.xticks([], [])  # nascondi read_name per chiarezza
    plt.ylabel("Number of edited bases")
    plt.xlabel("Reads (deduplicated)")
    plt.title(f"Distribution of editing types per read" + (f" - {status}" if status else ""))
    plt.legend(title="Edit type")
    plt.show()

    # --- 3. Stacked barplot frameshift vs non-frameshift ---
    # Raggruppa totale basi modificate per read
    df_shift = df_long.groupby(['read_name', 'frameshift']).agg({'count': 'sum'}).reset_index()

    # Pivot per stacked barplot
    df_pivot = df_shift.pivot(index='read_name', columns='frameshift', values='count').fillna(0)

    # Colori: frameshift True = red, False = green
    colors = {True: 'red', False: 'green'}
    df_pivot.plot(kind='bar', stacked=True, color=[colors[col] for col in df_pivot.columns], figsize=(8,5))
    plt.ylabel("Number of edited bases")
    plt.xlabel("Reads (deduplicated)")
    plt.title(f"Frameshift vs Non-frameshift" + (f" - {status}" if status else ""))
    plt.legend(title='Frameshift')
    plt.show()
