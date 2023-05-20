from typing import Set, List, Tuple
import functools
import time


def timer(func):
    """Decorate a function to determine the runtime of the function"""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        val = func(*args, **kwargs)
        end = time.perf_counter()
        work_time = end - start
        print(f'Runtime {func.__name__}: {round(work_time, 4)} seconds.')
        return val
    return wrapper


@timer
def minimizer_old(seq: str, w: int, k: int) -> Set[tuple]:
    """Generate a set of all (w,k)-minimizers with position index from a sequence"""

    minimizers = set()
    l = w + k - 1  # window length
    lseq = len(seq)  # sequence length

    for i in range(lseq - l + 1):  # loop of all windows
        # create a generator of all k-mers
        kmers = ((seq[j: j + k], j) for j in range(i, i + w))
        minimizers.add(min(kmers))  # find a minimizer and add in the set

    # find (u,k)-minimizer for all u < w at the beginning and the end of a sequence
    for i in range(k, l):  # loop of all windows, where u < w
        begin_kmers = ((seq[j: j + k], j) for j in range(i - k + 1))
        end_kmers = ((seq[j: j + k], j) for j in range(lseq - i, lseq - k + 1))
        minimizers.add(min(begin_kmers))
        minimizers.add(min(end_kmers))

    return minimizers


@timer
def generate_minimizers2(seq: str, w: int, k: int) -> Set[tuple]:
    """Generate set of all (w,k)-minimizers with position index from a sequence
     using array to save current list of k-mers"""

    minimizers = [(seq[:k], 0)]
    init_array = [(seq[:k], 0)]

    for i in range(1, w):  # find (u,k)-minimizer for all 1 < u < w at the beginning of a sequence
        init_array.append((seq[i: i + k], i))
        minimizers.append(min(init_array))

    for i in range(w, len(seq) - k + 1):  # find all (w,k)-minimizers
        del init_array[0]
        init_array.append((seq[i: i + k], i))
        minimizers.append(min(init_array))

    for i in range(w - 1):  # find (u,k)-minimizer for all u < w at the end of a sequence
        del init_array[0]
        minimizers.append(min(init_array))

    return set(minimizers)


@timer
def generate_minimizers_list(seq: str, w: int, k: int) -> List[Tuple]:
    """Generate set of all (w,k)-minimizers with position index from a sequence
     using array to save current list of k-mers"""

    init_array = [(seq[:k], 0)]

    minimizers = [(seq[:k], 0)]
    minimizers_set = {(seq[:k], 0)}

    # find (u,k)-minimizer for all 1 < u < w at the beginning of a sequence
    for i in range(1, w):
        m = (seq[i: i + k], i)
        init_array.append(m)
        minimizer = min(init_array)
        if minimizer not in minimizers_set:
            minimizers_set.add(minimizer)
            minimizers.append(minimizer)

    # find all (w,k)-minimizers
    for i in range(w, len(seq) - k + 1):
        del init_array[0]
        m = (seq[i: i + k], i)
        init_array.append(m)
        minimizer = min(init_array)
        if minimizer not in minimizers_set:
            minimizers_set.add(minimizer)
            minimizers.append(minimizer)

    # find (u,k)-minimizer for all u < w at the end of a sequence
    for i in range(w - 1):
        del init_array[0]
        minimizer = min(init_array)
        if minimizer not in minimizers_set:
            minimizers_set.add(minimizer)
            minimizers.append(minimizer)

    return minimizers

@timer
def generate_minimizers(seq: str, w: int, k: int) -> List[Tuple]:
    """Generate set of all (w,k)-minimizers with position index from a sequence
     using array to save current list of k-mers and additional variables to exclude repeating minimizers"""

    minimizers = [(seq[:k], 0)]
    init_array = [(seq[:k], 0)]
    prev_minimizer = (seq[:k], 0) # save previous minimizer

    for i in range(1, w):  # find (u,k)-minimizer for all 1 < u < w at the beginning of a sequence
        init_array.append((seq[i: i + k], i))
        current_minimizer = min(init_array) # calculate and save current minimizer
        if current_minimizer != prev_minimizer:
            minimizers.append(current_minimizer)
            prev_minimizer = current_minimizer

    for i in range(w, len(seq) - k + 1):  # find all (w,k)-minimizers
        del init_array[0]
        init_array.append((seq[i: i + k], i))
        current_minimizer = min(init_array)
        if current_minimizer != prev_minimizer:
            minimizers.append(current_minimizer)
            prev_minimizer = current_minimizer

    for i in range(w - 1):  # find (u,k)-minimizer for all u < w at the end of a sequence
        del init_array[0]
        current_minimizer = min(init_array)
        if current_minimizer != prev_minimizer:
            minimizers.append(current_minimizer)
            prev_minimizer = current_minimizer

    return minimizers
