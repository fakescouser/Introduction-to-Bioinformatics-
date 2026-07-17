from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Align import PairwiseAligner

# Ρυθμίσεις Entrez
Entrez.email = "up1084672@ac.upatras.gr"  # Βάλε το email σου εδώ
human_id = "P99999"

# --- ΒΗΜΑ 1: Fetching της ανθρώπινης πρωτεΐνης ---
print(f"Ανάκτηση αλληλουχίας για το ID: {human_id}...")
with Entrez.efetch(db="protein", id=human_id, rettype="fasta", retmode="text") as handle:
    human_record = SeqIO.read(handle, "fasta")
    human_seq = human_record.seq

print(f"Η αλληλουχία ανακτήθηκε: {human_record.description[:50]}...")

# --- ΒΗΜΑ 2: BLAST Query (Περιορισμός στα Mammalia) ---
print("Εκτέλεση BLASTP στα Mammalia (αυτό μπορεί να πάρει λίγη ώρα)...")
result_handle = NCBIWWW.qblast(
    program="blastp", 
    database="nr", 
    sequence=human_seq,
    entrez_query="Mammalia[Orgn]",
    hitlist_size=5  # Περιορίζουμε στα πρώτα 5 hits για ταχύτητα
)

# --- ΒΗΜΑ 3: Data Collection & Parsing ---
blast_record = NCBIXML.read(result_handle)
aligner = PairwiseAligner()

print("\n--- ΑΠΟΤΕΛΕΣΜΑΤΑ ΣΥΓΚΡΙΣΗΣ ---")
print(f"{'Accession':<15} | {'Description':<40} | {'Identity %'}")
print("-" * 75)

for alignment in blast_record.alignments:
    hit_id = alignment.accession
    hit_def = alignment.hit_def
    
    # Ανάκτηση της αλληλουχίας του "hit" για ακριβή σύγκριση
    with Entrez.efetch(db="protein", id=hit_id, rettype="fasta", retmode="text") as h_handle:
        hit_record = SeqIO.read(h_handle, "fasta")
        hit_seq = hit_record.seq

    # --- ΒΗΜΑ 4 & 5: Alignment & Identity Calculation ---
    # Κάνουμε pairwise alignment
    alignments = aligner.align(human_seq, hit_seq)
    best_alignment = alignments[0]
    
    # Υπολογισμός ταυτότητας βάσει του ορισμού της άσκησης:
    # Μετράμε πόσα αμινοξέα είναι ίδια στις αντίστοιχες θέσεις
    matches = 0
    # Χρησιμοποιούμε το μικρότερο μήκος για αποφυγή index error αν διαφέρουν
    min_len = min(len(human_seq), len(hit_seq))
    
    for i in range(min_len):
        if human_seq[i] == hit_seq[i]:
            matches += 1
            
    identity_percentage = (matches / len(human_seq)) * 100

    print(f"{hit_id:<15} | {hit_def[:40]:<40} | {identity_percentage:.2f}%")