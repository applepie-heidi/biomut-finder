from minimizers import timer, generate_minimizers
import sys
from SW import smith_waterman
from mapper import find_aligning_minimizers
import my_global


def reverse_dna(s):
    d = {'A': 'C', 'C': 'A', 'G': 'T', 'T': 'G'}
    res = ''.join([d[i] for i in s])
    return res[::-1]


@timer
def alignment(k=10, w=5, genr="ACTAGG", genc="GGGCGGCGGGTTTAAGGCGTTTCCTC"):
    mc = generate_minimizers(genc, w, k)
    mr = generate_minimizers(genr, w, k)
    mc_set = set([i[0] for i in mc])
    mr_set = set([i[0] for i in mr])
    # print('length of intersection:', len(mc_set & mr_set))

    common_minimizers = find_aligning_minimizers(mr, mc)
    # print(f'length of set of aligning minimizers: {len(common_minimizers)}')

    if len(common_minimizers) < 50:
        print("bad sequence")
        return 0

    ref = len(genc) * ['']
    prev_read = 100000000000000
    prev_gen = 10000000000000

    for i, minimizer in enumerate(common_minimizers):
        # print(f'length between two minimizers k = {minimizer[1] - prev_read}')
        # if minimizer[1] - (prev_read + k) > 0:
        # print(f'{i}: read, length between two minimizers k = {minimizer[1] - (prev_read + k)}')
        # pass
        if minimizer[2] - (prev_gen + k) > 0:
            # print(f'{i}: gen, length between two minimizers k = {minimizer[2] - (prev_gen + k)}')
            # if minimizer[1] - (prev_read + k) > 0 or minimizer[2] - (prev_gen + k) > 0:
            # if minimizer[1] - prev_v2 > (k + w) * 10:
            # print(f'stop on {i}th iteration because bad sequence')
            # return 0

            str_read = genr[prev_read + k: minimizer[1]]
            str_gen = genc[prev_gen + k: minimizer[2]]
            res_gen, res_read, start_gen, start_read = smith_waterman(str_gen,
                                                                      str_read)

            if res_gen and res_read:
                z = list(zip(res_gen, res_read))
                for ind, t in enumerate(z):
                    if t[0] == '-':
                        my_global.insertions = my_global.insertions + (
                        (ind + prev_gen + k + start_read, t[1]),)
                        del z[ind]
                res_gen_no, res_read_no = zip(*z)
                # res_gen_no = ''.join(res_gen_no)
                res_read_no = ''.join(res_read_no)
            else:
                continue

            # print local alignment
            '''
            print('lens_before:', 'gen:', len(str_gen), len(str_read), 'lens_after:', 'gen:', len(res_gen), len(res_read))
            print(i)
            print(res_gen)
            print(len(res_gen) * '|')
            print(res_read)
            print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
            print()
            '''
            # for ind, letter in enumerate(res_gen):

            for ind, letter in enumerate(res_read_no):
                ref[ind + prev_gen + k + start_read] = letter
        prev_read = minimizer[1] + k
        prev_gen = minimizer[2] + k

        for ind, letter in enumerate(minimizer[0]):
            ref[ind + minimizer[2]] = letter

    return ref
