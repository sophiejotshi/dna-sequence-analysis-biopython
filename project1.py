from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

def analyze_sequence(record):
    seq = str(record.seq).upper()
    length = len(seq)

    # Nucleotide counts
    a = seq.count("A")
    t = seq.count("T")
    g = seq.count("G")
    c = seq.count("C")
    n = seq.count("N")

    # Percentages
    at = (a + t) / length * 100
    gc = (g + c) / length * 100

    # Translation
    trimmed_length = length - (length % 3)
    trimmed_seq = seq[:trimmed_length]

    protein = Seq(trimmed_seq).translate(to_stop=True)
    protein_length = len(protein)


    # Stop codon analysis
    stop_codons = ["TAA", "TAG", "TGA"]
    stop_count = 0
    for i in range(0, length - 2, 3):
        codon = seq[i:i+3]
        if codon in stop_codons:
            stop_count += 1

    return {
        "sequence_id": record.id,
        "dna_length": length,
        "gc_content_percent": round(gc, 2),
        "at_content_percent": round(at, 2),
        "A_percent": round(a / length * 100, 2),
        "T_percent": round(t / length * 100, 2),
        "G_percent": round(g / length * 100, 2),
        "C_percent": round(c / length * 100, 2),
        "N_count": n,
        "protein_length": protein_length,
        "stop_codon_count": stop_count
    }

# ---------- Main pipeline ----------
input_fasta = "sars_cov_2.fasta"
results = []

for record in SeqIO.parse(input_fasta, "fasta"):
    results.append(analyze_sequence(record))

df = pd.DataFrame(results)

output_file = "dna_features_sars_cov_2.csv"
df.to_csv(output_file, index=False)

print("Analysis complete!")
print(df.head())
