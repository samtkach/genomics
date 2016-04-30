'''
written by tkach
burrows-wheeler transform for genome data
'''

from functools import partial
import sys

def bwt_radix(s):
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



def radix_key(text, value, step):
    return text[(value+step) % len(text)]

def old_count_matches_nerrors(bw, s, n):
    ranks,totals = rank(bw)
    fc = first_col(totals)

    if s[-1] not in fc:
        return None

    l,r = fc[s[-1]]
    possibles = range(l,r)
    errors = [0] * len(possibles)

    while i >= 0 and r > 1:
        char = s[i]

        for j in range(1,len(possibles)-1):
            ll = fc[char][0] + ranks[char][l-j]
            rr = fc[char][0] + ranks[char][r-j]
            



        l = fc[char][0] + ranks[char][l-1]
        r = fc[char][0] + ranks[char][r-1]
        i-=1
       
    return r-l

def rank_naive(bw):
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

def count_matches_naive(bw, s):
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

# -----------------------------------------------------------------------------




def suffix_array(s):
    '''build suffix array of string s'''
    sa = sorted([(s[i:], i) for i in xrange(0, len(s)+1)])
    return map(lambda x: x[1], sa)



def bwt(t):
    '''compute burrows-wheeler transform of string t'''
    bw = []
    for i in suffix_array(t):
        if i == 0:
            bw.append('$')
        else:
            bw.append(t[i-1])
    return ''.join(bw)




def first_col(totals):
    '''make dict of chars to range of occurences in first column'''
    col = {}
    temp = 0
    for i,j in sorted(totals.iteritems()):
        col[i] = (temp, temp+j)
        temp += j
    return col



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


def rank(bw):
    ''''''
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

#    print(totals) ##DEBUG
#    print(ranks) ##DEBUG
    return ranks, totals


def count_matches_exact(bw,s):
    '''return number of exact matches of s to bw'''
    ranks, totals = rank(bw)
    fc = first_col(totals)

    if s[-1] not in fc:
        return 0 #char isn't in bwt

    l,r = fc[s[-1]]
    i = len(s)-2
    while i >= 0 and r > 1:
        char = s[i]
        l = fc[char][0] + ranks[char][l-1] #R(aW)
        r = fc[char][0] + ranks[char][r-1] #Rbar(aW)
        i -= 1
        print('l: '+str(l)+' r: '+str(r))
        print(''.join(bw[l:r]))
        print(''.join(sorted(bw)[l:r]))

    return r-l


def inexact_search(bw, bwr, s, diff):
    '''find suffix array intervals with up to diff differences'''
    #ranks, totals 
    O,tot = rank(bw)

    #reverse ranks
    Oprime,junk = rank(bwr)
    
    #C[a] := number of lexicographically smaller letters than a in bw/reference
    C = compute_C(tot) 

    #D[i] := lower bound on number of differences in substring s[1:i] 
    D = compute_D(s, C, Oprime, bw)

    #call the secursive search function and return a list of SA-range tuples
    return inexact_recurse(s, len(s)-1, diff,1,len(bw)-1, D,C,O)

def compute_C(totals):
    '''compute C, the number of lexographicalls greater symbols in the ref'''
    C = {'a':0,'c':0,'g':0,'t':0,'$':0}
    for k in totals:
        for ref in totals:
            if ref < k: C[k] += totals[ref]

#    print(C) ##DEBUG
    return C


def compute_D(s, C, Oprime, bw):
    '''compute estimated lower bounds of differences in substring s[0:i] for all  in [0,len(s)]'''

    k = 1
    l = len(bw)-2
    z = 0
    D = [0] * len(s)

    for i in range(0, len(s)):
        k = C[s[i]] + Oprime[s[i]][k-1] + 1
        l = C[s[i]] + Oprime[s[i]][l]
        if k > l:
            k = 1
            l = len(bw)-1
            z = z+1
        D[i] = z
        print(str(k) + ', ' + str(l))

#    print('compute_D() returning:\n' + str(D)) ##DEBUG
    return D

def inexact_recurse(s,i,diff,k,l,D,C,O):
    '''recursion for inexact searching'''
    if diff < D[i]: return []
    
    if i < 0: return [(k,l)]

    I = []
    I = I + inexact_recurse(s,i-1,diff-1,k,l,D,C,O) #gap
    
    for b in ['a','c','g','t']:
        k = C[b] + O[b][k-1] + 1
        l = C[b] + O[b][l]

        if k<=l:
            I = I + inexact_recurse(s,i,diff-1,k,l,D,C,O) #gap
            if b == s[i]:
                I = I + inexact_recurse(s,i-1,diff,k,l,D,C,O) #match
            else:
                I = I + inexact_recurse(s,i-1,diff-1,k,l,D,C,O) #substitution

#    print('inexact_recurse returning:\n' + str(list(set(I)))) ##DEBUG
    return list(set(I))





def test():
    s = 'actactgacgtcagtcagttcagtacgtactactact'
    sa = suffix_array(s)
    bw = bwt(s)
    bwr = bwt(s[::-1])
    sa_ranges = inexact_search(bw,bwr,'act',0)
    for r in sa_ranges:
        (l,u) = r
        print('(l='+str(l)+',u='+str(u)+'):\n'+str(sa[l:u]))

        for i in sa[l:u]:
            print(s[i:]+', '+str(i))
    
    


 #   print(count_matches_exact(bw, 'atgatg'))
 #   print(count_matches_nerrors(bw,bwr,'atgatg',1)
    
test()
#print(bwt2(sys.argv[1]))
#print(suffix_array(sys.argv[1]))
