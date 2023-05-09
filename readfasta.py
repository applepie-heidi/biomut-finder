from Bio import SeqIO


def read_fasta(filename: str):
    record_list = []
    for record in SeqIO.parse(filename, "fasta"):
        record_list.append(record)
    return record_list
