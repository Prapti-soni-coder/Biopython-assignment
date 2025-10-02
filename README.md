# Biopython-assignment
#Prapti soni
# Install Biopython
 
!pip install --upgrade biopython

# -------------------------------
# Imports
# -------------------------------
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, seq1
from Bio.PDB import PDBList, PDBParser
from Bio.PDB.Polypeptide import PPBuilder, is_aa
from Bio.Align import PairwiseAligner
from Bio.Blast import NCBIWWW, NCBIXML
import os

# -------------------------------
# Settings
# -------------------------------
Entrez.email = "your.email@example.com"
NUC_ACC = "NM_001355006"
PDB_ID = "1TUP"
OUTDIR = "pipeline_output"
os.makedirs(OUTDIR, exist_ok=True)

# -------------------------------
# Fetch nucleotide sequence
# -------------------------------
handle = Entrez.efetch(db="nucleotide", id=NUC_ACC, rettype="fasta", retmode="text")
nuc_record = SeqIO.read(handle, "fasta")
handle.close()
nuc_seq = nuc_record.seq

# Translate full sequence (stop at first stop codon)
protein_from_nuc = nuc_seq.translate(to_stop=True)

# Shortened protein for alignment (first 200 aa)
protein_short = protein_from_nuc[:200]

print("\n=== Nucleotide Sequence ===")
print("ID:", nuc_record.id)
print("Length:", len(nuc_seq))
print("GC%: {:.2f}%".format(gc_fraction(nuc_seq) * 100))
print("Complement (first 50):", nuc_seq.complement()[:50])
print("Reverse Complement (first 50):", nuc_seq.reverse_complement()[:50])
print("Transcribed (first 50):", nuc_seq.transcribe()[:50])
print("Translated (first 50 aa):", protein_from_nuc[:50])

# -------------------------------
# Fetch protein structure (PDB)
# -------------------------------
pdbl = PDBList()
pdb_file = pdbl.retrieve_pdb_file(PDB_ID, pdir=OUTDIR, file_format="pdb")

parser = PDBParser(QUIET=True)
structure = parser.get_structure(PDB_ID, pdb_file)

# Extract first chain with standard amino acids
ppb = PPBuilder()
prot_seq = ""
for model in structure:
    for chain in model:
        seq_list = []
        for pp in ppb.build_peptides(chain):
            for res in pp:
                if is_aa(res, standard=True):
                    seq_list.append(seq1(res.get_resname()))
        if seq_list:
            prot_seq = "".join(seq_list)
            print(f"\n=== Protein Sequence from Chain {chain.id} ===")
            print("Length:", len(prot_seq))
            print("First 50 aa:", prot_seq[:50])
            prot_short = prot_seq[:200]  # Shortened for alignment
            break
    if prot_seq:
        break

# -------------------------------
# Pairwise Alignment
# -------------------------------
if len(protein_short) == 0 or len(prot_short) == 0:
    print("One of the sequences is empty. Alignment skipped.")
else:
    aligner = PairwiseAligner()

    # Global alignment
    aligner.mode = "global"
    global_alignment = aligner.align(protein_short, prot_short)
    if len(global_alignment) > 0:
        print("\n=== Global Alignment ===")
        print(global_alignment[0])
    else:
        print("No global alignment found.")

    # Local alignment
    aligner.mode = "local"
    local_alignment = aligner.align(protein_short, prot_short)
    if len(local_alignment) > 0:
        print("\n=== Local Alignment ===")
        print(local_alignment[0])
    else:
        print("No local alignment found.")

# -------------------------------
# Find first start codon
# -------------------------------
for i in range(0, len(nuc_seq)-2, 3):
    if nuc_seq[i:i+3] == "ATG":
        print("\nFirst start codon found at nucleotide position:", i)
        break

# -------------------------------
# Perform BLAST (protein)
# -------------------------------
print("\nSubmitting BLAST (blastp) request... this may take several minutes")
blast_file = os.path.join(OUTDIR, "blast_results.xml")
result_handle = NCBIWWW.qblast("blastp", "nr", str(protein_from_nuc), hitlist_size=5)

with open(blast_file, "w") as f:
    f.write(result_handle.read())
result_handle.close()

# Parse BLAST results
with open(blast_file) as f:
    blast_records = NCBIXML.read(f)

print("\n=== Top 5 BLAST Hits ===")
if len(blast_records.alignments) == 0:
    print("No BLAST hits found.")
else:
    for i, alignment in enumerate(blast_records.alignments[:5]):
        print(f"\nHit {i+1}:")
        print("Title:", alignment.title)
        print("Length:", alignment.length)
        print("E-value:", alignment.hsps[0].expect)

print("\n=== Pipeline Finished ===")

