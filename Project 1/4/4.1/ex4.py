from Bio import SeqIO
import os

# --- CONFIGURATION ---
# Using the paths from your desktop directory
input_fasta = r"C:\Users\User\Desktop\bioinfo\3\SARS-CoV-2_sequence.fasta"
output_fasta = r"C:\Users\User\Desktop\bioinfo\3\SARS-CoV-2_spike_sequence.fasta"

def extract_spike_protein(input_path, output_path):
    print(f"Searching for Spike protein (gene=S) in: {input_path}")
    
    found = False
    try:
        with open(input_path, "r") as infile:
            # Parse the main fasta file
            for record in SeqIO.parse(infile, "fasta"):
                # Identify the spike protein by the "gene=S" tag in the header
                if "gene=S" in record.description:
                    with open(output_path, "w") as outfile:
                        SeqIO.write(record, outfile, "fasta")
                    
                    found = True
                    print("--- Extraction Complete ---")
                    print(f"ID: {record.id}")
                    print(f"Description: {record.description}")
                    print(f"Sequence Length: {len(record.seq)} amino acids")
                    print(f"Saved to: {output_path}")
                    break
        
        if not found:
            print("Status: Process finished, but no entry with 'gene=S' was found.")
            
    except FileNotFoundError:
        print(f"Error: The file '{input_path}' was not found. Please check the path.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    extract_spike_protein(input_fasta, output_fasta)