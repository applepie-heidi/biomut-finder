import sys

from minimizers import generate_minimizers
from readfasta import read_fasta
from mapper import find_aligning_minimizers


def reversed_complement(sequence: str) -> str:
    """Return reversed complement of a sequence"""
    complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join(complement[base] for base in reversed(sequence))


def main_test():
    gen = "AGGGGGAGTTTGAGA"
    gen2 = "AGCTAGTGTTTGAGAATTTTTTGTACCGCGTACGTTGCACCGTACCCAGTCTCTCGCGCGTGCGTCAACGTACGTCGAGACTGCATGCATGCGCGTGCAGTTTTTTTTCGCGCGCGCGCGCGCGCGCGTGGGTGTGTGTGTGTGTGTGTGTGTGCGCGCGCGCGCGCGCAATATTTAAATTTCCCGGGCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTCGCGCGCGCGCGGGGGAAGGGAGGGGAGGGAGGGAGGAGGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAAAAAAAACCCCTCTTTCTTTTCACACCCACCCCACCCCACCCACCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCTCTTTTTTTTTTTTACCCACACACACACACACACACACAGGGTGTGTGTGTGTGTGTGTGTGTGTGCGCGTCGCTCTCTCTCTCTCCCCTTTCCCGGGAAATTCCGCTGCGTCGTCGTCGTCGTCGTAAACGCGCGCGCGCCCGGGTTTGGCCCGGTTCCAAGCGCGCGTGTGTGTGAACGCGCGCGCGTGTGGAGAGAGAGAGAGTTTTTTTTTTTTGTGTACCCCGCGTTGTCTCTCTGCGAACGATCGCGCGCGCGCGCGCGCTGGGGTGGGGGTGGGGTGGGGCGGGGCGCGAAAAAAAAAAAAAGAGAGGGGGGGGTGTGTGTGAGTGACCCGTCAGTCAGTCGTACGTTGGGGGGGTGTGGCAGTGGGGGGGGGGGGGGGGGGGCGTGTGTGTGGCGTGCAGTCAGTCAGTCAGTCAGTACGTACGTACGTACCCCCAGTGAGAGTGGGGTCAGTCAAAAGTCAGTTTTGCGTTGGGTACGTACGTACGTTGCATGGGCGGGGGGAAAAACCCCCTTTTGGGCGTGGTTTTTGTGGGGTGTGTGTGTGTGTGGGGGGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCGCGCGCTGCGTCGTCGTCGTCTGCTGCTGCCGTCGAAATGGCCAAATGCATGCCCACGTCATGCATGCATGCATGCATGCATGCGTACGTACGTCATGCGTGGGCGTACGTCGTGCGTGGCGTGGCCCAATGC"
    k, w = 4, 6
    mg1_1 = generate_minimizers(gen, w, k)
    mg1_2 = generate_minimizers(gen, w, k)
    print(mg1_1)
    print(mg1_2)
    print(mg1_1 == mg1_2)
    print(
        f'memory for minimizer: {sys.getsizeof(mg1_1)}, memory for minimizer_second: {sys.getsizeof(mg1_2)}: ')

    mg2_1 = generate_minimizers(gen2, w, k)
    mg2_2 = generate_minimizers(gen2, w, k)
    print(mg2_1)
    print(mg2_2)
    print(mg2_1 == mg2_2)
    print(
        f'memory for minimizer: {sys.getsizeof(mg2_1)}, memory for minimizer_second: {sys.getsizeof(mg2_2)}: ')


def main():
    k, w = 10, 5

    gen_ref = read_fasta("data/lambda.fasta")[0].seq
    gen_reads = read_fasta("data/lambda_simulated_reads.fasta")

    # turn into string IMPORTANT, 7 times faster
    gen_ref = str(gen_ref)

    ref_minimizers = generate_minimizers(gen_ref, w, k)

    for gen_read in gen_reads:
        gen_read = str(gen_read.seq)
        gen_read_reversed = reversed_complement(gen_read)
        read_minimizers = generate_minimizers(gen_read, w, k)
        read_minimizers_reversed = generate_minimizers(gen_read_reversed, w, k)
        print(find_aligning_minimizers(ref_minimizers, read_minimizers))
        print(find_aligning_minimizers(ref_minimizers, read_minimizers_reversed))
        print()


if __name__ == '__main__':
    main()
