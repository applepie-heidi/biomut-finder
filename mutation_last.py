from Bio import SeqIO
from collections import Counter
from minimizers import timer
import alignment_last
import my_global
import csv


def count_letters(ref_letter='T',
                  tuple_letters=('A', 'A', 'T', 'A', '-', 'T', 'A')):
    c = Counter(tuple_letters)
    max_letter = (max(c.items(), key=lambda x: x[1]))
    if max_letter[0] == ref_letter:
        return 0
    else:
        max_count = max_letter[1]
        all_count = len(tuple_letters)
        pr = max_count / all_count
        if max_letter[0] == '-' or all_count > 1:
            # if all_count > 1:
            if pr >= 0.7:
                # print(f'mutation: {max_letter[0]} instead {ref_letter}, {c}, {tuple_letters}')
                # return 1
                # return ref_letter
                return max_letter[0]
            else:
                return 0


def reverse_dna(s):
    d = {'A': 'C', 'C': 'A', 'G': 'T', 'T': 'G'}
    res = ''.join([d[i] for i in s])
    return res[::-1]


@timer
def find_mutations(k=10, w=5):
    data = []
    seq = []
    for record in SeqIO.parse("lambda.fasta", "fasta"):
        seq.append(record.seq)

    res_seq = [tuple(i, ) for i in seq[0]]

    reads = []
    for record in SeqIO.parse("lambda_simulated_reads.fasta", "fasta"):
        reads.append(record.seq)

    print("alingment:")
    genc = str(seq[0])
    bads = []
    for i, genr in enumerate(reads):
        # if i > 200:
        # break
        bds = [17, 18, 21, 22, 24, 26, 27, 31, 33, 34, 35, 36, 37, 38, 40, 41,
               46, 47, 49, 50, 51, 52, 55, 58, 59, 60, 63, 64, 66, 67, 68, 73,
               74, 75, 77, 79, 82, 83, 84, 88, 89, 92, 93, 94, 95, 96, 97, 98,
               100, 102, 103, 104, 105, 106, 107, 108, 110, 111, 116, 117, 119,
               120, 122, 123, 124, 125, 127, 128, 129, 130, 131, 132, 133, 134,
               135, 136, 137, 138, 140, 141, 142, 145, 147, 149, 151, 152, 153,
               154, 155, 156, 160, 161, 165, 166, 168, 169, 171, 172, 175, 177,
               178, 180, 181, 184, 185, 186, 188, 189, 190, 197, 199, 201, 202,
               203, 211, 213, 214, 216, 220, 222, 223, 224, 226, 228, 229, 231,
               232, 235, 236, 237, 241, 242, 244, 245]
        if i < 15 or i in bds:
            # if 15 >= i or i >= 17:
            continue

        genr = str(genr)
        print(f'#{i}', len(genc), len(genr))
        # count += alignment.alignment(k, w, genr, genc)
        res = alignment_last.alignment(k, w, genr, genc)
        # ins = Counter(my_global.insertions)
        # print('ins', len(my_global.insertions), ins)
        # print("---------------------reversed read-------------------------")
        # count += alignment.alignment(k, w, reverse_dna(genr), genc)
        # print(res)

        if res:
            for ind, letter in enumerate(res):
                if letter != '':
                    res_seq[ind] = res_seq[ind] + (letter,)
        else:
            bads.append(i)
        print("-----------------------------------------------------")

    empty_number = 0
    ins = Counter(my_global.insertions)
    for i, j in ins.items():
        if j > 2:
            # res_seq[i[0]] = res_seq[i[0]] + (('I', i[0], i[1]),) * j
            # print(res_seq[i[0]])
            # print(f'mutation insertion: I, {i[0]}, {i[1]}')
            data.append(['I', i[0], i[1]])

    for ind, i in enumerate(res_seq):
        if len(i) != 1:
            res = count_letters(i[0], i[1:])
            if res:
                # print(i, ind)
                # print('---------------------')
                # mutation = ['X', ind, res]
                if res == '-':
                    mutation = ['D', ind, res]
                else:
                    mutation = ['X', ind, res]
                # print(mutation)
                data.append(mutation)
        else:
            empty_number += 1

    # print(empty_number)
    # print(data)

    columns = ['mutation', 'position', 'nucleotide']
    with open('mutations.csv', 'w', encoding='utf-8', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(columns)
        for row in sorted(data, key=lambda x: x[1]):
            writer.writerow(row)
