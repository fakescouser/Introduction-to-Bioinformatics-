from Bio import AlignIO
import pandas as pd
import os

def analyze_alignment(filepath, method_name, format="fasta"):
    if not os.path.exists(filepath):
        print(f"  Warning: File '{filepath}' not found. Skipping...")
        return None

    try:
        alignment = AlignIO.read(filepath, format)
    except Exception as e:
        print(f" Error reading {filepath}: {e}")
        return None

    n_seqs = len(alignment)
    aln_len = alignment.get_alignment_length()

    total_positions = n_seqs * aln_len
    total_pairs = n_seqs * (n_seqs - 1) / 2

    identity_matches = 0
    total_comparisons = 0
    gap_count = 0
    sp_score = 0

    for i in range(aln_len):
        column = alignment[:, i]
        for a in range(n_seqs):
            for b in range(a + 1, n_seqs):
                aa = column[a]
                bb = column[b]

                if aa == "-" or bb == "-":
                    gap_count += 1
                    continue
                total_comparisons += 1

                if aa == bb:
                    identity_matches += 1
                    sp_score += 1

    # Calculate metrics
    percent_identity = (identity_matches / total_comparisons) * 100 if total_comparisons > 0 else 0
    # Normalized by all possible sequence pairs across all positions
    percent_gaps = (gap_count / (total_positions * (n_seqs - 1))) * 100 
    normalized_sp = (sp_score / (total_pairs * aln_len))  

    # CHANGED: All rounding updated to 4 decimal places
    return {
        "Method": method_name,
        "Total Seqs": n_seqs,
        "Length": aln_len,
        "Percent Identity": round(percent_identity, 4),
        "Percent Gaps": round(percent_gaps, 4),
        "Normalized SP Score": round(normalized_sp, 4)
    }

# --- THE EXECUTION PART ---
if __name__ == "__main__":
    files_to_run = [
        ("mafft_defaultnew2.aln-fasta", "MAFFT Default"),
        ("mafft_custom1.aln-fasta", "MAFFT Custom"),
        ("mafft_custom_gappenalty3_ext1.aln-fasta", "MAFFT Custom 2"),
        ("muscle_default.aln-fasta", "MUSCLE Default"),
        ("muscle_first.aln-fasta", "MUSCLE First"),
        ("muscle_second.aln-fasta", "MUSCLE Second"),
        ("tcoffee_default.aln-fasta", "T-COFFEE Default"),
        ("tcoffee_blosum_alligned.aln-fasta", "T-COFFEE Blosum Alligned"),
        ("tcoffee_paminput.aln-fasta", "T-COFFEE PAM Input")
    ]
    
    all_results = []

    for file_path, label in files_to_run:
        res = analyze_alignment(file_path, label)
        if res:
            all_results.append(res)

    if all_results:
        df = pd.DataFrame(all_results)
        
        print("\n" + "="*70)
        print("             FINAL ALIGNMENT COMPARISON (High Precision)")
        print("="*70)
        # Using to_string ensures the console doesn't truncate the decimals
        print(df.to_string(index=False))
        print("="*70)

        try:
            df.to_excel("Alignment_Comparison_Detailed.xlsx", index=False)
            print("\n High-precision results saved to 'Alignment_Comparison_Detailed.xlsx'!")
        except Exception as e:
            print(f"\nCould not save to Excel: {e}")
    else:
        print("No files were processed. Check your filenames!")