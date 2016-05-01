'''
sam tkach, kenny xiao
inexact search over burrows-wheeler genome data
'''

import sys
from enum import Enum
from operator import itemgetter

from bwt import *

# bw is bwt of genome (B)
# bwr is bwt of reverse genome (B')
# s is short read to be matched (W)
# diff is max num of differences (z)

#alphabet of symbols in allowed
alphabet = set(['A', 'T', 'C', 'G'])

C = {} # 
O = {} # defined below
D = [] #

# enum represents the type of the alignment choice
Type = Enum('Type', 'START MATCH MISMATCH INSERTION DELETION')

# rewards/penalties:
gap_open = 0
gap_ext = 1
mismatch = 1
match = 0

# option switches:
NO_INDELS = False

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

    #call the recursive search function and return a list of SA-range tuples
    sa_index_set = inexact_recursion(s, len(s)-1, diff,0,len(bw)-1, Type.START)
    index_dict = {}

    for (i,j) in sa_index_set:
        #if index already exists, pick the higher diff value
        if i in index_dict:
            if index_dict[i] < j:
                index_dict[i] = j
        else:
            index_dict[i] = j

    # sort list by diff from highest to lowest
    return sorted(index_dict.items(), key=itemgetter(1), reverse=True) 



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

def inexact_recursion(s,i,diff,k,l, prev_type):
    '''search bwt recursively and tolerate errors'''

    #pruning based on estimated mistakes
    if diff < get_D(i):
        return set()

    #end of query condition
    temp = set()
    if i < 0:
        for j in range(k,l+1):
            temp.add((j,diff))
        return temp

    #search
    sa_idx = set() #set of suffix array indices at which a match starts 
    
    if not NO_INDELS:
    # Insertion
        if prev_type == Type.INSERTION:
            sa_idx = sa_idx.union(inexact_recursion(s,i-1,diff-gap_ext,k,l,Type.INSERTION))
        else:
            sa_idx = sa_idx.union(inexact_recursion(s,i-1,diff-gap_ext-gap_open,k,l,Type.INSERTION))

    for char in alphabet:
        temp_k = C[char] + get_O(char, k-1) + 1
        temp_l = C[char] + get_O(char, l)

        if temp_k <= temp_l:
            if not NO_INDELS:
            # Deletion
                if prev_type == Type.DELETION:
                    sa_idx = sa_idx.union(inexact_recursion(s,i,diff-gap_ext,temp_k,temp_l,Type.DELETION))
                else:
                    sa_idx = sa_idx.union(inexact_recursion(s,i,diff-gap_ext-gap_open,temp_k,temp_l,Type.DELETION))
            if char == s[i]:
                # Match!
                sa_idx = sa_idx.union(inexact_recursion(s,i-1,diff+match,temp_k,temp_l,Type.MATCH))
                
            else:
                # Mismatch
                sa_idx = sa_idx.union(inexact_recursion(s,i-1,diff-mismatch,temp_k,temp_l, Type.MISMATCH))

    #print diff
    return sa_idx


def print_output(sa_index_list, sa, s):
    '''print formatted output'''
    sa_values = [(sa[i],j) for (i,j) in sa_index_list]


    print '-----------------------------------------'
    print str(len(sa_values)) + " match(es) found!\n"
    print "Score\t\tPosition\tSuffix\n"
    for v,x in sa_values:
        print str(x) + "\t\t" + str(v) + "\t\t" + s[v:]

    print '----------------------------------------'

def test():
    s = 'ATGCGTAATGCCGTCGATCG'
    sa = suffix_array(s)
    bw = bwt(s)
    bwr = bwt(s[::-1])

    ##DEBUG
    #print("BW: " + bw) 
    #print("BWR: " + bwr + '\n')

    print_output(inexact_search(bw,bwr,'GTA',1), sa, s)
    #print "\nfinal ranges: " + str(sa_ranges) + "\n"


def main():

    threshold = 1 # this will be the z value
    usage = ('\nusage: python search_bwt.py [--no-indels] [test|<reference file name>] [<read file name>]\n')

    if '--no-indels' in sys.argv:
        global NO_INDELS
        NO_INDELS = True


    if len(sys.argv) == 1:
        print usage
        return

    elif sys.argv[1].lower() == 'test':
        test()
        return
    

    elif len(sys.argv) < 3:
        print usage
        return

    fread = open(sys.argv[-1])
    fref = open(sys.argv[-2])

    ref = ''.join(fref.readlines()).replace('\n','')
    read = ''.join(fread.readlines()).replace('\n','')

    sa = suffix_array(ref)

    bw = bwt(ref)
    bwr = bwt(ref[::-1])

    print read
    print_output(inexact_search(bw,bwr,read,threshold), sa, ref)

    fread.close()
    fref.close()


main()
