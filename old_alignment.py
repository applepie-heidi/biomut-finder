from minimizers import minimizer, minimizer_second, timer, generate_minimizers_new
from NW import nw
from LIS import find_lis, find_aligning_minimizers

def reverse_dna(s):
    d = {'A': 'C', 'C': 'A', 'G': 'T', 'T': 'G'}
    res = ''.join([d[i] for i in s])
    return res[::-1]


@timer
def alignment(k=10, w=5, genr="ACTAGG", genc="GGGCGGCGGGTTTAAGGCGTTTCCTCTTCGTCATAACT"):
    mc = generate_minimizers_new(genc, w, k)
    mr = generate_minimizers_new(genr, w, k)
    mc_set = set([i[0] for i in mc])
    mr_set = set([i[0] for i in mr])
    # print('intersection:', mc_set & mr_set)
    print('length of intersection:', len(mc_set & mr_set))

    common_minimizers = find_aligning_minimizers(mr, mc)
    #print(f'result of aligning : {common_minimizers}')
    print(f'length of set of aligning minimizers: {len(common_minimizers)}')
    ##print("--------------------------------------------------------------")

    if len(common_minimizers) < 20:
        print("bad sequence")
        return 0

    ref = len(genc) * [''] # empty sequence
    prev = 100000000000000
    prev_r = 10000000000000

    # find dynamic alingment for gaps
    for i, minimizer in enumerate(common_minimizers):
        if minimizer[1] - prev > k:
            ##if minimizer[1] - prev > (k + w) * 10:
            ##    print(f'stop on {i}th iteration because bad sequence')
            ##    return 0

            str_r = genr[prev + k: minimizer[1]]
            str_g = genc[prev_r + k: minimizer[2]]
            res_g, res_r = nw(str_g, str_r) #global dynamic algorithm

            print(i)
            print(res_g)
            print(len(res_g) * '|')
            print(res_r)
            print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

            #for ind, letter in enumerate(res_g): #!!!!!!!!!!!!
            for ind, letter in enumerate(res_r):
                # ref[ind + prev_r_v2 + k] = letter
                # print(ind, letter)
                    if ref[ind + prev_r + k] == '':
                        ref[ind + prev_r + k] = letter + '@' #'@' for check, that this is dynamic alingment
                        #ref[ind + prev_r + k] = letter
                    else:
                        #ref[ind + prev_r + k] = ref[ind + prev_r + k] + '&' + letter + '@'
                        ref[ind + prev_r + k] = ref[ind + prev_r + k] + letter + '&' # '&' check that we have rewriting (no good)
            #TODO add for the second alignment

        prev = minimizer[1] # add k ?
        prev_r = minimizer[2] # add k ?

        # minimizers
        for ind, letter in enumerate(minimizer[0]):
            # print(ind, letter)
            if ref[ind + minimizer[2]] == '':
                ref[ind + minimizer[2]] = letter
            else:
                ref[ind + minimizer[2]] = ref[ind + minimizer[2]] + '#' + letter # check rewriting using '#'
            # ref[ind + minimizer[2]] = letter
            # ref[ind + minimizer[2]] = (letter, ind + minimizer[1], ind + minimizer[2], i)
    #print(ref) #check result
    return ref
