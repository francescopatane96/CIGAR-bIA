# 🧬 CIGAR-bIA  
**CIGAR-based Indels Analysis**

<p align="center">
  <img src="cigar_bia/logo.png" alt="CIGAR-bIA Logo" width="350" height="350/>
</p>



CIGAR-bIA is a Python toolkit for analyzing genomic editing events based on CIGAR strings extracted from BAM files.  
It provides functions to identify, quantify, and visualize insertions and deletions across specific genomic regions.

---

## ⚙️ Installation

### 🔹 From terminal
```bash
git clone https://github.com/francescopatane96/CIGAR-bIA.git
cd CIGAR-bIA
pip install .
```

### 🔹 From Jupyter Notebook
```python
!git clone https://github.com/francescopatane96/CIGAR-bIA.git
!pip install CIGAR-bIA/.
```

---

## 💻 Example Usage

### 1️⃣ Import the main modules
```python
import cigar_bia.events as cba
import cigar_bia.plotting as cbp
```

### 2️⃣ Analyze indel events from CIGAR strings
```python
events = cba.analyze_editing_events(
    bam_file="/data/possorted_genome_bam.bam",
    chrom="chr20",
    start=40688387,
    end=40688661,
    meta_file="/results/meta_mafb_KO.csv",
    status="MAFB-KO"
)
```

### 3️⃣ Visualize the results
```python
cbp.plot_cigar_reads(events, chrom="chr20", start=40688387, end=40688661)
```

---

## 📄 License
This project is released under the **MIT License**.  
See the [LICENSE](./LICENSE) file for more details.

---

## 🤝 Contributing
Contributions, bug reports, and feature suggestions are welcome!  
Please open a **pull request** or submit a **GitHub issue** on the project page.

---

📘 *Author:* [Francesco Patanè](https://github.com/francescopatane96)
