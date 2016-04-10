'''
written by tkach
burrows-wheeler transform for genome data
'''

from functools import partial
import sys

def bwt2(s):
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



def bw2_key(text, value, step):
    return text[(value+step) % len(text)]



def suffix_array(s):
    sa = sorted([(s[i:], i) for i in xrange(0, len(s)+1)])
    return map(lambda x: x[1], sa)



def bwt(t):
    bw = []
    for i in suffix_array(t):
        if i == 0:
            bw.append('$')
        else:
            bw.append(t[i-1])
    return ''.join(bw)



def rank(bw):
    '''counts total occurences of each character and 
       counts number of occurences of a row's char until that row (rank)'''
    totals = dict()
    ranks = []
    for char in bw:
        if char not in totals.keys():
            totals[char] = 0
        ranks.append(totals[char])
        totals[char] += 1
    return ranks, totals

def first_col(totals):
    '''make dict of chars to range of occurences in first_col'''
    first = {}
    temp = 0
    for i,j in sorted(totals.iteritems()):
        first[i] = (temp, temp+j)
        temp += j
    return first



def ibwt(bw):
    '''decode bwt'''
    ranks, totals = rank(bw)
    fc = first_col(totals)
    row = 0
    t = "$"

    while bw[row] != '$':
        char = bw[row]
        t = char+t
        row = fc[char][0] + ranks[row]

    return t



def count_matches(bw, s):
    '''see how many matches s has in bw'''
    ranks,totals = rank(bw)
    fc = first_col(totals)
    l,r = fc[s[-1]]
    i = len(s)-2
    while i >= 0 and r > 1:
        char = s[i]
        # iterate and find occurences of char
        j = 1
        while j < r:
            if bw[j] == char: 
                l = first[char][0] + ranks[j]
                break
            j+=1

        # no occurences condition
        if j == r: 
            l = r
            break

        r-=1
        while bw[r] != char:
            r-=1

        r = fc[char][0] + ranks[r] + 1
        i-=1

    return r - 1



def rank_O1(bw):

    totals = {}
    ranks = {}

    for char in bw:
        if char not in totals:
            totals[char] = 0
            ranks[char] = []

    for char in bw:
        totals[char] += 1
        for t in totals.iterkeys():
            ranks[t].append(totals[t])

    return ranks, totals



def count_matches_O1(bw,s):
    ranks, totals = rank_O1(bw)
    fc = first_col(totals)

    if s[-1] not in fc:
        return 0 #char isn't in bwt

    l,r = fc[s[-1]]
    i = len(s)-2
    while i >= 0 and r > 1:
        char = s[i]
        l = fc[char][0] + ranks[char][l-1]
        r = fc[char][0] + ranks[char][r-1]
        i -= 1
    return r-l




def test():
    bw = bwt(sys.argv[1])
    print(ibwt(bw))
    print(bw)
    print(count_matches_O1(bw, sys.argv[2]))

test()
#print(bwt2(sys.argv[1]))
#print(suffix_array(sys.argv[1]))
