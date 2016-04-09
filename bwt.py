'''
written by tkach
burrows-wheeler transform for genome data
'''

from functools import partial
import sys

def bwt(s):
    s += "\0"
    return ''.join(s[i-1] for i in radix_sort(range(len(s)), partial(bw_key, s)))



def radix_sort(values, key, step=0):
    if len(values) < 2:
        for v in values:
            yield v
        return

    bins = {}
    for v in values:
        bins.setdefault(key(v,step), []).append(v)

    for k in sorted(bins.keys()):
        for r in radix_sort(bins[k], key, step + 1):
            yield r


def bw_key(text, value, step):
    return text[(value+step) % len(text)]


print(bwt(open(sys.argv[1]).readlines()))
