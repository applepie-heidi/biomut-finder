from collections.abc import Sequence
from typing import Tuple, List


# import math


def find_lis(x):  # Sequence[Tuple]) -> List[Tuple]:
    """Find the longest increasing subsequence.
Description of the algorithm and pseudo-code at the link:
https://en.wikipedia.org/wiki/Longest_increasing_subsequence#Efficient_algorithms"""

    n = len(x)
    p = [0] * n
    m = [0] * (n + 1)
    # m[0] = -1
    l = 0

    for i in range(n):
        left = 1
        right = l
        while left <= right:
            # mid = math.ceil((left + right)/2)
            mid = (left + right) // 2
            if x[m[mid]][2] < x[i][2]:
                left = mid + 1
            else:
                right = mid - 1
        newl = left
        p[i] = m[newl - 1]
        m[newl] = i

        if newl > l:
            l = newl

    s = [0] * l
    k = m[l]
    for i in range(l - 1, -1, -1):
        s[i] = x[k]
        k = p[k]
    return s


def find_aligning_minimizers(first_minimizers, second_minimizers):
    # Sequence[Tuple], second_minimizers: Sequence[Tuple]) -> List[Tuple]:
    """Find list of  minimizers for aligning with index in first sequence 
       and index in second sequence using LIS algorithm"""

    res = []
    for i1, j1 in first_minimizers:
        for i2, j2 in second_minimizers:
            if i1 == i2:
                res.append((i1, j1, j2))
    return find_lis(res)


def score_alignment(sequence1, sequence2):
    """Score alignment of two sequences"""

    score = 0
    for i in range(len(sequence1)):
        if sequence1[i] == sequence2[i]:
            score += 1
    return score


def align_sequences(sequence1: str, sequence2: str, aligning_minimizers):
    """Align two sequences using aligning minimizers"""

    subsequences = []
    for i, j1, j2 in aligning_minimizers:
        # extend alignment
        start_position = j1 - j2
        end_position = j1 + len(sequence2) - j2
        subsequence = sequence1[start_position: end_position]
        score = score_alignment(subsequence, sequence2)
        subsequences.append((subsequence, start_position, score))
    # find best alignment
    return max(subsequences, key=lambda x: x[2])
