import collections

from mapper import find_aligning_minimizers, align_sequences
from minimizers import generate_minimizers
from readfasta import read_fasta


def reversed_complement(sequence: str) -> str:
    """Return reversed complement of a sequence"""
    complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join(complement[base] for base in reversed(sequence))


def find_mutations(gen_ref: str, sequence, filename: str):
    """Find mutations in a sequence compared to a reference genome
    and write them to a file"""
    with open("data/bla.txt", "w") as f:
        for s in sequence:
            f.write(s + "\n")

    with open(filename, "w") as file:
        for i in range(len(gen_ref)):
            if sequence[i]:
                sequence_list = sequence[i].split(",")
                chosen, chosen_count = \
                    collections.Counter(sequence_list).most_common(1)[0]
                if chosen != gen_ref[i]:
                    if chosen == "-":
                        file.write(f"D,{i},-\n")
                    elif len(chosen) > 1:
                        for letter in range(chosen):
                            file.write(f"I,{i},{letter}\n")
                    else:
                        file.write(f"X,{i},{chosen}\n")
                    file.write(f"{i} {chosen} {chosen_count}\n")


def test_main():
    k, w = 10, 5
    gen_ref = read_fasta("data/lambda.fasta")[0].seq
    gen_read = "GGGACGAAGGTTCCTGCGCGGTTAGATCGATTATTATTCGGTGCCGATATTCGTAGAACAAAACCGAGGACACTGGAAGGGCGAGTAGCTCTCTCTCGCAAGGGGACTTGGGTACAAGGATGCAAACCTGCAAAAACCAATATCAGGCTGGGGCTACATTTTCGAGTAGGAATACACAGTCGTGTAACCTTGAGATCTTTCACATTAATCCCAACGACCACAGTCGTACACGAGCATCTTACGGAAGTTTATTTAGATTCCGGGCAAAGCGCATCAACAGCGGGTATGTTGTAATAAACCGCTTGGGGTGATTTTGGAGCCATAAAGCTGGCACTTGTTGCGAGTCGCTTTGGGGTGAAGTCTCATACCTCTCGGGTCGGGACTGCCAGACAGGACTAGCAGCCCTCGTCGCCAATATCATCGCCCCAAGGCAGTAACATATTCCAATAGTTAGTGGCGCTTAATTCGATACCGTAACAGATTACTACTACAATTCGGTATGTCATGATGACTGATCGTCACCTGGATAAAGTGTGTTATCTAGGCGCGTCTCTGAGGTCAGCTATTTCCTTTTTTGCTTAAGAAACCTCCGGGAACCAGTGCACCGCCGGACGTGAAGTGCAGGAGCACAATCCGCACTGCCCGACACGCTCAAGCAGATGAGCGGCGCTCCCCGGTTTCACGGATATGTGGTTGGTGGCGACTATTAGTATTAGACCCTTCGTCGTACAGTTTATCCGGAGCGTAGCGGGGCGGTTCGCGATGCTGACTTCCTTCGCCATGTGCCCGTGCGGGATAGAAAGGTCTAATGACATATGGTTTAAGACCTATACGGATCGAGGCACTTAGACCAAGACGCCTTGTCTCAAATCTACCCATAGCAACACCCTGGACCATCGTGAATAACCTCTTAGTGACCAATACTGATCAGCGGCGATAGTCGCCGGGCCTGCAAGCTGACACGTGGACGCATCGCGGAACCGTTTTCAGGGCCCTACAGGCAAAAAATAGAGTCGGGCGTAAGCCGGCGATCGCAATACCTCTATAAACGAGTACGACCTGCTGACTCCTGTGGGCCCCAACCTTTTCCACATTGGCCAATTGTATTTACGGTATTCTTATACAAACTTGCTGGTTGGCTCTGTGACCGTCCTGAGCATTCAATTACAGCCACGAAGCACTACTATTGCCATCCTGCACACCCCGATAGATAGCCTATCAAGAGACCTTGAATCGCCCTACCAAAGCACTTTGCAAGTACTCTGACCTAGGGAATCATAAGTCTAGCGCTGCGATATACCCATTACCTCCGGCGGAGGAGGCGGTGCGATCTAAAAAATGTAGTATAAACATCAAAGATAGTCAAGGATGTTGTTTGACACACGGGGTCAGATTAGGATTAGTTCTCGGAGAGATGCCTTGTCCAACTCAGACTAAAGGGCAAGCCATACAGTGATATAGCTACCCCGGGGGCTACCAAGTCGTTCATGTCCGAGGGGGTTTCACAATCCATAGGGACCTGATCGACCCTTTCCCACTCTAGCCGAGCTTTTGGGATCTTGTGCCGTTACGGATGACGCGTGGGCTTTGATGATAATTGCGAGGTGGGCGTCATCCCGTCAGACCAGGACGAATTTTCACCAACAAAATGGGCGTCGTATCTTTGGCATTTGTGGAATGGACCAACCGGTAATGGGGCGATCAACTAGACTGCTACCGCGCCATGCATGATATTAGAAGTCCGGGGGGCTGGCAAGCCATTAGTGATTTGCAGCGTGGGGGACGACCGCGGATAAGCCGCCTCTCCTTATGGGTCCGAGGCAATCTCCCCAGCTATAGTCCTGAATCCCATCTTTAATCACTGACGAAGCGGTCGCGTTACCTCGCCGAGCTGCCAGCTACTATGTCACGCGTCACCCTGCCTCCGTGCATTGGAAAAAACGCTCAGAACAAGTACTCTGTTACGCAACACCATCCAAGAAGGCAATTTGCTAGCAATCGGACCGATGTGGGTGGTAGTGAGCACGTCTGTGCCCCTCGGGATTAAATGATACGGCCTTTACGGTGAGACACGGAAGGTCGCTTTCACCAGTTAGTGAAGCTACTTTCGCCACGCTGCGAAGGATATTTTGCCATACCGCACTAACTCCGTTGCCACCTAGGTTCACAAAAGATAATGAGGCATTACCCCGACTACATTTACGGGGCAAGGATTATAAAGCTATGGGCTCGTCAATAGTCCCAAACGATCCAAGGTGTTTGGCAACTAAGCATCATATACCACGAGATCGCCAGCGCTTCCCCGGAGAATAAATAAAGTGCTTTAAGTGGTCACTCCGTTGAAGTAAGTCCGCATATAGCTGCAGGGGGGGGTAGCTTCCCTGTACGGGCCCCCATATCTACAGGTGCATATAACTTCCGCGACCGGGTCGCTACTATTATAGTCCGCAAAGCCCAGAACGGGTCTCCGTGTGAATATGTGGGAGCAGCCGGAAAACATTAAAAAACGGAGTGTATCATCATAGGTGACTAAGCGCCCGAGAGACATCAAGTGGTAAAGACACTAATCCGTGGAAGTATGTACTGCCCGGAGGAGACTACGCATGGTAGAACCATACAATGCACCTAAGTTCTCAGATCGCGAACCTTATAATTATATAACTCGGTACAGTTACGCGTTAAATGAGACAATTCGACATACTGAGGTTGCGGACCTATGTGCACGCTTGACTGCTGTCGACTAGCTGTGGTGCTGCCTGTGCGGGTGGGATTATGAAGCCATGGCGGTGCGGGCTCGAAGTGTCTCGATGGCCCCCGCTCACTTTCGACCACGCAGATGATATAATCATCATTTCTCAGCTTGATGGTGTCAAATCCTAAGAGAGAAACTGGGGCCAAGATGCGAGAGAATACACGGGCAGTCGGCAAAGAGAATTCTCAGTGTCTGTTCAAGGATGTTGTCTCATTGCCAACAGACCTGAGGACTCCTCGCCTTCAAATTTTTTGCGCTCTGCCCCCCGAAAACTTGAGTGAGTAGGTCGGG"
    ref_minimizers = generate_minimizers(str(gen_ref), w, k)
    sequence = [""] * len(gen_ref)
    read_minimizers = generate_minimizers(gen_read, w, k)
    aligning_minimizers = find_aligning_minimizers(ref_minimizers,
                                                   read_minimizers)
    start_position, result, score = align_sequences(str(gen_ref), gen_read,
                                                    aligning_minimizers)
    print(start_position, result, score)


def main():
    k, w = 10, 5
    mini_n = 20
    gen_ref = read_fasta("data/ecoli.fasta")[0].seq
    gen_reads = read_fasta("data/ecoli_simulated_reads.fasta")

    gen_ref = read_fasta("data/lambda.fasta")[0].seq
    gen_reads = read_fasta("data/lambda_simulated_reads.fasta")

    # turn into string IMPORTANT, 7 times faster
    gen_ref = str(gen_ref)

    ref_minimizers = generate_minimizers(gen_ref, w, k)

    sequence = [""] * len(gen_ref)
    for gen_read in gen_reads:
        gen_read = str(gen_read.seq)
        gen_read_reversed = reversed_complement(gen_read)
        read_minimizers = generate_minimizers(gen_read, w, k)
        read_minimizers_reversed = generate_minimizers(gen_read_reversed, w, k)

        start_position, result, score = 0, [], 0
        start_position_reversed, result_reversed, score_reversed = 0, [], 0
        if len(read_minimizers) >= mini_n:
            aligning_minimizers = find_aligning_minimizers(ref_minimizers,
                                                           read_minimizers)
            start_position, result, score = align_sequences(gen_ref, gen_read,
                                                            aligning_minimizers)
        if len(read_minimizers_reversed) >= mini_n:
            aligning_minimizers_reversed = find_aligning_minimizers(
                ref_minimizers,
                read_minimizers_reversed)
            start_position_reversed, result_reversed, score_reversed = align_sequences(
                gen_ref, gen_read_reversed,
                aligning_minimizers_reversed)

        if score != 0 or score_reversed != 0:
            try:
                if score >= score_reversed:
                    for i in range(len(result)):
                        sequence[start_position + i] += result[i] + ","
                else:
                    for i in range(len(result_reversed)):
                        sequence[start_position_reversed + i] += result_reversed[
                                                                     i] + ","
            except IndexError:
                raise IndexError("Index out of range")
    find_mutations(gen_ref, sequence, "data/mutations.txt")


if __name__ == '__main__':
    main()
