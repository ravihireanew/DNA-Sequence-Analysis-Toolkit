# 🧬 DNA Sequence Analysis Toolkit

A powerful, interactive bioinformatics web application for analyzing DNA sequences. Built with **Biopython** and **Streamlit**.

This tool allows users to fetch sequences from NCBI, perform comprehensive analysis, visualize results, and export data — all through an easy-to-use web interface.



## ✨ Features

### Core Analysis
- Fetch DNA sequences directly from **NCBI** using accession numbers
- Support for **local FASTA file upload** and raw sequence pasting
- Calculate **GC Content**, Molecular Weight, and Nucleotide Composition
- Find **Open Reading Frames (ORFs)** in all 6 reading frames
- Generate **Reverse Complement**

### Advanced Visualizations
- Interactive **GC Content Sliding Window** plot
- **ORF Length Distribution** histogram
- **Nucleotide Composition** pie chart
- **Codon Usage** analysis (Top 20 codons)

### Additional Tools
- Restriction enzyme site finder (EcoRI, HindIII, BamHI, etc.)
- Self Dot Plot for sequence repeats
- Safe handling of ambiguous bases (`N`)
- Summary statistics and downloadable results

## 🚀 Live Demo
(Coming soon — will be added after deployment on Streamlit Cloud)

## 📸 Screenshots

*(Add 2–3 screenshots here: Overview tab, ORF tab, and GC plot)*

## 🛠️ Technologies Used

- **Python 3.10+**
- **Streamlit** — Web interface
- **Biopython** — Sequence analysis
- **Plotly** — Interactive visualizations
- **Pandas** — Data handling

## 📥 Installation & Setup

### 1. Clone the Repository
git clone https://github.com/ravihireanew/dna-sequence-analysis-toolkit.git
cd dna-sequence-analysis-toolkit

📊 Project Highlights

Handles both short mRNA and long genomic sequences
Robust error handling for ambiguous bases (N)
Professional interactive UI
Ready for portfolio & academic projects

🚀 Future Enhancements (Planned)

Protein translation & amino acid analysis
Multiple sequence alignment & comparison
Primer design tool
BLAST integration
Codon Adaptation Index (CAI)
One-click deployment on Streamlit Cloud
