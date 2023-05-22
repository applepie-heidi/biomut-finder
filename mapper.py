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


def smith_waterman_algorithm(sequence1, sequence2):
    match_score = 2
    mismatch_penalty = -1
    gap_penalty = -2

    # Initialize the scoring matrix
    rows = len(sequence1) + 1
    cols = len(sequence2) + 1
    scoring_matrix = [[0] * cols for _ in range(rows)]

    # Initialize variables to store the maximum score and its position
    max_score = 0
    max_position = (0, 0)

    # Fill the scoring matrix
    for i in range(1, rows):
        for j in range(1, cols):
            # Calculate the match/mismatch score
            if sequence1[i - 1] == sequence2[j - 1]:
                match = scoring_matrix[i - 1][j - 1] + match_score
            else:
                match = scoring_matrix[i - 1][j - 1] + mismatch_penalty

            # Calculate the gap scores
            delete = scoring_matrix[i - 1][j] + gap_penalty
            insert = scoring_matrix[i][j - 1] + gap_penalty

            # Determine the maximum score for the current position
            score = max(0, match, delete, insert)
            scoring_matrix[i][j] = score

            # Update the maximum score and its position if necessary
            if score > max_score:
                max_score = score
                max_position = (i, j)

    # Traceback to find the local alignment
    alignment = []
    i, j = max_position
    while i > 0 and j > 0 and scoring_matrix[i][j] > 0:
        if scoring_matrix[i][j] == scoring_matrix[i - 1][j - 1] + (
                match_score if sequence1[i - 1] == sequence2[
                    j - 1] else mismatch_penalty):
            alignment.append((sequence1[i - 1], sequence2[j - 1]))
            i -= 1
            j -= 1
        elif scoring_matrix[i][j] == scoring_matrix[i - 1][j] + gap_penalty:
            alignment.append((sequence1[i - 1], '-'))
            i -= 1
        else:
            alignment.append(('-', sequence2[j - 1]))
            j -= 1

    alignment.reverse()  # Reverse the alignment list to obtain the correct order

    return alignment, max_score


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
