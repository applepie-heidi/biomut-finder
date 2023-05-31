import alignment
from Bio import SeqIO
from collections import Counter
import alignment_new

def count_letters(ref_letter='T', tuple_letters=('A', 'A', 'T', 'A', '-', 'T', 'A')):
    """ determine the most frequent letter and is a mutation or not """

    c = Counter(tuple_letters)
    max_letter = (max(c.items(), key=lambda x: x[1]))
    if max_letter[0] == ref_letter:
        #print('no mutations')
        return 0
    else:
        max_count = max_letter[1]
        all_count = len(tuple_letters)
        pr = max_count / all_count
        if pr >= 0.7:
            print(f'mutation: {max_letter[0]} instead {ref_letter}, {c}, {tuple_letters}')
            return 1
        else:
            #(print('no mutation anyway'))
            return 0

def reverse_dna(s):
    d = {'A': 'C', 'C': 'A', 'G': 'T', 'T': 'G'}
    res = ''.join([d[i] for i in s])
    return res[::-1]

def find_mutation():

    k, w = 10, 5
    count = 0
    seq = []  # Setup an empty list
    for record in SeqIO.parse("lambda.fasta", "fasta"):
        seq.append(record.seq)
    res_seq = [tuple(i,) for i in seq[0]]
    #res_seq = [tuple()] * len(seq[0])

    reads = []  # Setup an empty list
    for record in SeqIO.parse("lambda_simulated_reads.fasta", "fasta"):
        reads.append(record.seq)

    print("alingment:")

    genc = str(seq[0]) # reference sequence
    bads_seqs = []
    for i, genr in enumerate(reads): # loop for all reads
        # if i > 200: # exclude "bad" reads to accelerate code
           # break
        # bds = [17, 18, 21, 22, 24, 26, 27, 31, 33, 34, 35, 36, 37, 38, 40, 41, 46, 47, 49, 50, 51, 52, 55, 58, 59, 60, 63, 64, 66, 67, 68, 73, 74, 75, 77, 79, 82, 83, 84, 88, 89, 92, 93, 94, 95, 96, 97, 98, 100, 102, 103, 104, 105, 106, 107, 108, 110, 111, 116, 117, 119, 120, 122, 123, 124, 125, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 140, 141, 142, 145, 147, 149, 151, 152, 153, 154, 155, 156, 160, 161, 165, 166, 168, 169, 171, 172, 175, 177, 178, 180, 181, 184, 185, 186, 188, 189, 190, 197, 199, 201, 202, 203, 211, 213, 214, 216, 220, 222, 223, 224, 226, 228, 229, 231, 232, 235, 236, 237, 241, 242, 244, 245]
        # if i < 15 or i in bds:
            # continue

        genr = str(genr) #read
        print(f'#{i}', len(genc), len(genr))
        #count += alignment.alignment(k, w, genr, genc) # old version don't work
        res = alignment_new.alignment(k, w, genr, genc)
        #print("---------------------reversed read-------------------------")
        #count += alignment.alignment(k, w, reverse_dna(genr), genc)
        #print(res)

        if res:
            for ind, letter in enumerate(res):
                if letter != '':
                    res_seq[ind] = res_seq[ind] + (letter,) #aling our reed to reference sequence (actually add in tuple new letter)
                #print(letter)
        else:
            bads_seqs.append(i) # find number of "bad" sequence for future

    # print(res_seq)

    # find mutations
    '''
    for ind, i in enumerate(res_seq):
        if len(i) != 1:
            res = count_letters(i[0], i[1:])
            if res:
                print(i, ind)
                print('---------------------')
    '''