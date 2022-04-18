import Bio
from Bio import SeqIO

def readFasta(fasta_file):
    """
    Convert fasta file to list of tuples
    """
    fasta_list = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        fasta_list.append([record.id, record.description, record.seq.lower()])
    return fasta_list

def write_fasta(fasta_record, fasta_file):
    """
    Write list of tuples to fasta file
    """
    with open(fasta_file, 'w') as f:
        f.write('>' + fasta_record[0] + ' ' + fasta_record[1] + '\n')
        f.write(str(fasta_record[2]) + '\n')