'''
written by tkach
burrows-wheeler transform for genome data
'''

import sys

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
    '''rank(char) := list of number of occurences of a char for each substring R[:i] (reference)'''
    totals = {}
    ranks = {}


    for char in alphabet:
        if (char not in totals) and (char != '$'):
            totals[char] = 0
            ranks[char] = []


    for char in bw:
        if char != '$':
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


# bw is bwt of genome (B)
# bwr is bwt of reverse genome (B')
# s is short read to be matched (W)
# diff is max num of differences (z)

alphabet = set(['A', 'T', 'C', 'G'])
C = {}
O = {}
D = []

def inexact_search(bw, bwr, s, diff):
    '''find suffix array intervals with up to diff differences'''
    
    
    #ranks, totals
    #O is a dictionary with keys $,A,C,G,T, and values are arrays of counts
    global O
    O,tot = rank(bw)
    
    ##DEBUG
    #print('This is O:\n')
    #print O

    #reverse ranks
    Oprime,junk = rank(bwr)
    
    ##DEBUG
    #print('This is Oprime:\n')
    #print Oprime

    #C[a] := number of lexicographically smaller letters than a in bw/reference
    global C
    C = compute_C(tot) 

    #D[i] := lower bound on number of differences in substring s[1:i] 
    global D
    D = compute_D(s, C, Oprime, bw)

    #call the secursive search function and return a list of SA-range tuples
    return inexact_recursion(s, len(s)-1, diff,0,len(bw)-1)


def compute_C(totals):
    '''compute C, the number of lexographically greater symbols in the ref'''
    C = {'A':0,'C':0,'G':0,'T':0}
    for k in alphabet:
        for ref in alphabet:
            if ref < k: C[k] += totals[ref]

#    print('compute_C() returning: ' + str(C) + '\n') ##DEBUG
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
        #print(str(k) + ', ' + str(l))

#    print('compute_D() returning: ' + str(D) + '\n') ##DEBUG
    return D



def get_D(i):
    '''enforse condiion that if D[i] is set to -1, its 
       value will be considered as 0'''
    if i < 0:
        return 0
    else:
        return D[i]

def get_O(char, index):
    '''see get_D()'''
    if index < 0:
        return 0
    else:
        return O[char][index]

def inexact_recursion(s,i,diff,k,l):
    '''search bwt recursively and tolerate errors'''

    #pruning based on estimated mistakes
    if diff < get_D(i):
        return set()

    #end of query condition
    temp = set()
    if i < 0:
        for j in range(k,l+1):
            temp.add(j)
        return temp

    #search
    sa_idx = set() #sset of suffix array indicces at which a match starts 
    sa_idx = sa_idx.union(inexact_recursion(s,i-1,diff-1,k,l))
    for char in alphabet:
        temp_k = C[char] + get_O(char, k-1) + 1
        temp_l = C[char] + get_O(char, l)

        if temp_k <= temp_l:
            sa_idx = sa_idx.union(inexact_recursion(s,i,diff-1,temp_k,temp_l))
            if char == s[i]:
                sa_idx = sa_idx.union(inexact_recursion(s,i-1,diff,temp_k,temp_l))
            else:
                sa_idx = sa_idx.union(inexact_recursion(s,i-1,diff-1,temp_k,temp_l))

    return sa_idx


def test():
    s = 'ATGCGTAATGCCGTCGATCG'
    sa = suffix_array(s)
    bw = bwt(s)
    bwr = bwt(s[::-1])

    ##DEBUG
    #print("BW: " + bw) 
    #print("BWR: " + bwr + '\n')

    sa_index_set = inexact_search(bw,bwr,'GTA',1)
    sa_values = [sa[i] for i in sa_index_set]

    print str(len(sa_values)) + " match(es) found!\n"
    print "Position\tSuffix\n"
    for v in sa_values:
        print str(v) + "\t\t" + s[(v-1):]

    #print "\nfinal ranges: " + str(sa_ranges) + "\n"


def main():
    

    threshold = 100 # this will be the z value

    if len(sys.argv) == 1:
        print ('usage: python bwt.py [test|<reference file name>] [<read file name>]')
        return

    elif sys.argv[1].lower() == 'test':
        test()
        return
    

    elif len(sys.argv) < 3:
        print ('usage: python bwt.py [test|<reference file name>] [<read file name>]')
        return

    fread = open(sys.argv[2])

    ref = fref.readlines().replace('\n','')
    read = fread.readlines().replace('\n','')

    bw = bwt(ref)
    bwr = bwt(ref[::-1])

    sa = suffix_array(ref)

    sa_indices = inexact_search(bw,bwr,read,threshold)

    sa_values = [sa[i] for i in sa_indices]

    print str(len(sa_values)) + " match(es) found!\n"
    print "Position\tSuffix\n"
    for v in sa_values:
        print str(v) + "\t\t" + s[(v-1):]



main()
