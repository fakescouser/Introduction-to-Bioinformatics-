from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# Βήμα 1: Ορισμός της αλληλουχίας
sequence = "GIVEQCCTSICSLYQLENYCN"

print("Εκτέλεση BLAST... Παρακαλώ περιμένετε (μπορεί να πάρει 5-6 λεπτά).")

# Βήμα 2: Χρήση της qblast 
# Αν θέλουμε αναζήτηση ΜΟΝΟ για ανθρώπους, βάζουμε # στην αρχή των 2ον και 3ον φίλτρων.
# Αν θέλουμε τη βασική αναζήτηση, βάζουμε ένα # στην αρχή των 1ον και 3ον.
# Αν θέλουμε να περιορίσουμε την αναζήτηση μόνο σε "Bacteria" χρησιμοποιούμε το 3ο φίλτρο και βάζουμε # στα άλλα δύο παραπάνω

#result_handle = NCBIWWW.qblast("blastp", "nr", sequence, entrez_query="Homo sapiens[ORGN]") # Φίλτρο Ανθρώπου
#result_handle = NCBIWWW.qblast("blastp", "nr", sequence) # Βασική αναζήτηση
result_handle = NCBIWWW.qblast("blastp", "nr", sequence, entrez_query="Bacteria[ORGN]") # Bacteria 

# Βήμα 3: Ανάγνωση αποτελεσμάτων
blast_record = NCBIXML.read(result_handle)

# Εύρεση του πρώτου hit (με το χαμηλότερο E-value)
if blast_record.alignments:
    first_alignment = blast_record.alignments[0]
    hsp = first_alignment.hsps[0]
    
    print("\n--- Αποτελέσματα BLAST ---")
    print(f"Πρωτεΐνη/Τίτλος: {first_alignment.title}")
    print(f"E-value: {hsp.expect}")
else:
    print("Δεν βρέθηκαν αποτελέσματα.")