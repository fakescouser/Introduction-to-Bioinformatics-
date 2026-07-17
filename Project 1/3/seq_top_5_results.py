from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# Βήμα 1: Ορισμός της αλληλουχίας
sequence = "GIVEQCCTSICSLYQLENYCN"

print("Εκτέλεση BLAST... Παρακαλώ περιμένετε.")

# Βήμα 2: Εκτέλεση qblast
# Μπορείς να ενεργοποιήσεις/απενεργοποιήσεις το φίλτρο ανθρώπου με το '#'
#result_handle = NCBIWWW.qblast("blastp", "nr", sequence, entrez_query="Homo sapiens[ORGN]")
result_handle = NCBIWWW.qblast("blastp", "nr", sequence)

# Βήμα 3: Ανάγνωση αποτελεσμάτων
blast_record = NCBIXML.read(result_handle)

# Ορίζουμε πόσα hits θέλουμε να εμφανίσουμε
top_hits = 5 

print(f"\n--- Εμφάνιση των πρώτων {top_hits} αποτελεσμάτων ---")

# Χρησιμοποιούμε slicing [:top_hits] για να πάρουμε τα πρώτα 5 alignments
for i, alignment in enumerate(blast_record.alignments[:top_hits]):
    hsp = alignment.hsps[0] # Παίρνουμε το πρώτο High-scoring Segment Pair
    print(f"\nHit #{i + 1}:")
    print(f"Πρωτεΐνη/Τίτλος: {alignment.title}")
    print(f"E-value: {hsp.expect}")
    print(f"Score: {hsp.score}")
    print("-" * 30)

if not blast_record.alignments:
    print("Δεν βρέθηκαν αποτελέσματα.")