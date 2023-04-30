from collections.abc import Sequence
# import math


def find_lis(x: Sequence[int]) -> list[int]:
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
            if x[m[mid]] < x[i]:
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
