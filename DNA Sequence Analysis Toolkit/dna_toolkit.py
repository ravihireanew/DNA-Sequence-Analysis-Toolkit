from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, molecular_weight
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from datetime import datetime

# Always set your email when using NCBI (required!)
Entrez.email = "kumbharavindra0@gmail.com"   # ← CHANGE THIS TO YOUR EMAIL

def fetch_sequence(accession: str, db="nucleotide"):
    """Fetch DNA sequence from NCBI using accession number (e.g., NC_3479)"""
    print(f"Fetching sequence {accession} from NCBI...")
    try:
        handle = Entrez.efetch(db=db, id=accession, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        print(f"Success! Sequence length: {len(record.seq)} bp")
        return record
    except Exception as e:
        print(f"Error fetching sequence: {e}")
        return None

def analyze_sequence(record):
    """Perform all basic analyses on a SeqRecord"""
    seq = record.seq
    dna_seq = seq.upper() if not seq.isupper() else seq

    results = {
        "Accession": record.id,
        "Description": record.description,
        "Length_bp": len(dna_seq),
        "GC_Content_%": round(gc_fraction(dna_seq) * 100, 2),
        "Molecular_Weight_DNA_g_mol": round(molecular_weight(dna_seq.replace("N", ""), seq_type="DNA"), 2) if "N" in dna_seq else round(molecular_weight(dna_seq, seq_type="DNA"), 2),
    }
    # Reverse Complement
    rev_comp = dna_seq.reverse_complement()

    # Translate in 3 forward frames (simple)
    protein_frames = []
    for frame in range(3):
        length = 3 * ((len(dna_seq) - frame) // 3)
        protein = dna_seq[frame:frame + length].translate(to_stop=False)
        protein_frames.append(str(protein))

    # Find simple ORFs (longest ORF in all 6 frames)
    orfs = find_orfs(dna_seq)
    results["Longest_ORF_aa"] = max(len(orf) for orf in orfs) if orfs else 0
    results["Number_of_ORFs"] = len(orfs)

    print("\n=== Analysis Results ===")
    for key, value in results.items():
        print(f"{key:20}: {value}")

    return {
        "seq": dna_seq,
        "rev_comp": rev_comp,
        "proteins": protein_frames,
        "orfs": orfs,
        "results": results
    }

def find_orfs(dna_seq, min_protein_length=50, table=11):
    """Find potential Open Reading Frames in all 6 frames"""
    orfs = []
    for strand, nuc in [(+1, dna_seq), (-1, dna_seq.reverse_complement())]:
        for frame in range(3):
            length = 3 * ((len(nuc) - frame) // 3)
            pro = nuc[frame:frame + length].translate(table)
            pro_str = str(pro)
            # Split by stop codon (*) and keep long ones
            for p in pro_str.split("*"):
                if len(p) >= min_protein_length:
                    orfs.append(p)
    return orfs

def visualize_results(record, analysis):
    """Create nice visualizations"""
    seq = analysis["seq"]
    
    # GC content sliding window
    window = 100
    gc_values = [gc_fraction(seq[i:i+window]) * 100 for i in range(0, len(seq)-window, window//2)]
    
    plt.figure(figsize=(12, 8))
    
    plt.subplot(2, 2, 1)
    plt.plot(gc_values)
    plt.title("GC Content Sliding Window (100 bp)")
    plt.xlabel("Position")
    plt.ylabel("GC %")
    
    plt.subplot(2, 2, 2)
    lengths = [len(orf) for orf in analysis["orfs"]]
    if lengths:
        sns.histplot(lengths, bins=20, kde=True)
        plt.title("ORF Length Distribution")
        plt.xlabel("Length (amino acids)")
    
    plt.subplot(2, 2, 3)
    plt.text(0.1, 0.5, f"Sequence Length: {len(seq)} bp\n"
                       f"GC Content: {analysis['results']['GC_Content_%']}%\n"
                       f"Longest ORF: {analysis['results']['Longest_ORF_aa']} aa",
             fontsize=12, bbox=dict(facecolor='lightblue', alpha=0.5))
    plt.axis('off')
    
    plt.tight_layout()
    plt.savefig(f"dna_analysis_{record.id}_{datetime.now().strftime('%Y%m%d')}.png", dpi=300)
    plt.show()

# ====================== MAIN ======================
if __name__ == "__main__":
    # Example accessions (change as you like)
    accessions = ["NM_000207.3", "AH002844.2"]   # pPCP1 plasmid and an example viral genome
    
    for acc in accessions:
        record = fetch_sequence(acc)
        if record:
            analysis = analyze_sequence(record)
            visualize_results(record, analysis)
            print("-" * 60)