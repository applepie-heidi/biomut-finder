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


def smith_waterman_algorithm(sequence1, sequence2):
    """Smith-Waterman algorithm for local alignment of two sequences"""

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

    # Reverse the alignment list to obtain the correct order
    alignment.reverse()

    return alignment, max_score


def decide_indels(subsequence1: str, subsequence2: str):
    score = 0
    if len(subsequence1) > len(subsequence2):
        result = ["-"] * len(subsequence1)
        i = -1
        for j, letter in enumerate(subsequence2):
            try:
                i = subsequence1[i + 1:].index(letter)
                if len(subsequence1) - i <= len(subsequence2) - j:
                    result[i] = letter
                    score += 1
                else:
                    result[j] = letter
                    i = j
            except ValueError:
                result[j] = letter
                i = j
    elif len(subsequence1) < len(subsequence2):
        result = []
        i = -1
        for j, letter in enumerate(subsequence1):
            try:
                i2 = subsequence2[i + 1:].index(letter)
                if len(subsequence2) - i2 <= len(subsequence1) - j:
                    result.append(subsequence2[i + 1:i2 + 1])
                    i = i2
                else:
                    result.append(letter)
                    i = j
                    score += 1
            except ValueError:
                result.append(letter)
                i = j
                score += 1
    else:
        result = []
        for i in range(len(subsequence1)):
            letter1 = subsequence1[i]
            letter2 = subsequence2[i]
            result.append(letter2)
            if letter1 == letter2:
                score += 1
    return result, score


def align_sequences(sequence1: str, sequence2: str, aligning_minimizers):
    """Align two sequences using aligning minimizers"""

    if len(aligning_minimizers) == 0:
        return -1, [], 0

    final_result = []
    final_score = 0

    i_prev, position1_prev, position2_prev = aligning_minimizers[0]
    start_position = position1_prev - position2_prev
    if start_position < 0:
        start_position = 0
    subsequence1 = sequence1[start_position:position1_prev]
    subsequence2 = sequence2[:position2_prev]
    result, score = decide_indels(subsequence1, subsequence2)
    final_result.extend(result)
    final_score += score

    for index in range(1, len(aligning_minimizers)):
        i, position1, position2 = aligning_minimizers[index]
        subsequence1 = sequence1[position1_prev:position1]
        subsequence2 = sequence2[position2_prev:position2]
        result, score = decide_indels(subsequence1, subsequence2)
        final_result.extend(result)
        final_score += score

        if index == len(aligning_minimizers) - 1:
            end_position = position1 + len(sequence2) - position2
            if end_position > len(sequence1):
                end_position = len(sequence1)
            subsequence1 = sequence1[position1:end_position]
            subsequence2 = sequence2[position2:]
            result, score = decide_indels(subsequence1, subsequence2)
            final_result.extend(result)
            final_score += score
        i_prev, position1_prev, position2_prev = i, position1, position2

    return start_position, final_result, final_score
