def minimizer(seq: str, w: int, k: int) -> set:

    '''Generate set of all (w,k)-minimizers with index of position from a sequence'''

    minimizers = set() 
    l = w + k - 1 # lenght of window
    lseq = len(seq) # lenght of sequence

    for i in range(lseq - l + 1): # loop of all windows
        kmers = ((seq[j: j + k], j) for j in range(i, i + w)) # create generator of all k-mers
        minimizers.add(min(kmers)) # find minimizer and add in set

    # find (u,k)-minimizer for all u < w at the beginning and end of a sequence
    for i in range(k, l): # loop of all windows, where u < w
        begin_kmers = ((seq[j: j + k], j) for j in range(i - k + 1))
        end_kmers = ((seq[j: j + k], j) for j in range(lseq - i, lseq - k + 1))
        minimizers.add(min(begin_kmers))
        minimizers.add(min(end_kmers))

    return minimizers