import sys
import random
from search_bwt import *
import time

def parse_reads(myfile):
    # Get tuples of identifier, position, and read

    file_name = myfile
    read_file = open(file_name, 'r')
    #new_file_name = file_name[:-6] + '_new.fastq'
    #write_file = open(new_file_name, 'w')

    line_num = 0

    read_dict = {}
    identifier = ''
    position = 0
    read = ''      

    for line in read_file:

        # increment line number
        line_num += 1

        # Get 1/100 of the reads (127)
        if line_num % 202 in [1,2]:

            #Get name and position
            if line_num % 2 == 1:

                word_arr = line.split()
                identifier = word_arr[0]
                position = int((word_arr[3])[4:])

            if line_num % 2 == 0:
                base = ''
                rand = random.randint(1,4)
                if rand == 1:
                    base = 'A'
                elif rand == 2:
                    base = 'T'
                elif rand == 3:
                    base = 'G'
                else:
                    base = 'C'

                read_dict[identifier] = (position, line.rstrip().replace('N',base))
                #write_file.write(identifier + "\t" + str(position) + "\n" + line.rstrip() + "\n")

    read_file.close()
    #write_file.close()

    return read_dict

def reverse_complement(seq):
    seq_dict = {'A':'T','T':'A','G':'C','C':'G'}
    return "".join([seq_dict[base] for base in reversed(seq)])  

def align_reads(genome_file, reads_file):

    # Get reads
    read_dict = parse_reads(reads_file)

    """
    for identifier in read_dict:
        print identifier
        print "pos=" + str(read_dict[identifier][0])
        print read_dict[identifier][1]
        print "----------------------------------------------------"
    """

    # Get reference genome

    fref = open(genome_file)
    ref = ''.join(fref.readlines()).replace('\n','')

    print "Calculating suffix array..."

    sa = suffix_array(ref)
    print "Calculating bw..."
    bw = bwt(ref)
    print "Calculating bwr..."
    bwr = bwt(ref[::-1])

    # Search over all reads

    print "Processing reads..."
    # key = identifier, value = best match position
    new_read_dict = {}
    threshold = 4

    
    count = 0
    total = len(read_dict)

    print "Progress (threshold=" + str(threshold) + "): "
    start = time.time()

    for identifier in read_dict:
        count += 1
        if count % 5 == 0: print "{0:.1f}".format(100.0*count/total) + "%"

        s = read_dict[identifier][1]
        best_position1, score1 = best_match_position(bw, bwr, s, threshold, sa)
        best_position2, score2 = best_match_position(bw, bwr, reverse_complement(s), threshold, sa)

        if score1 >= score2 and score1 != -1:
            new_read_dict[identifier] = best_position1
        elif score1 < score2:
            new_read_dict[identifier] = best_position2
        elif score1 == score2 and score1 == -1:
            new_read_dict[identifier] = -1              
        else:
            new_read_dict[identifier] = -2

    end = time.time()
    print "Done aligning reads (" + str(end-start) + " seconds elapsed)..."

    num_correct = 0
    num_incorrect = 0
    num_missing = 0
    num_error = 0
    for identifier in read_dict:
        if read_dict[identifier][0] == new_read_dict[identifier]:
            num_correct += 1
        elif new_read_dict[identifier] == -1:
            num_missing += 1
        elif new_read_dict[identifier] == -2:
            num_error += 1
        elif read_dict[identifier][0] != new_read_dict[identifier]:
            num_incorrect += 1

    print "Number of reads: \t\t" + str(len(read_dict))
    print "Number of correct alignments: \t" + str(num_correct)
    print "Number of incorrect alignments: " + str(num_incorrect)
    print "Number of 'no matches': \t" + str(num_missing)
    print "Number of errors: \t\t" + str(num_error)



    """
    if score1 >= score2 and score1 != -1:
        print "identifier: " + identifier + "\t\test. best position: " + \
        str(best_position1) + "\t\tactual position: " + str(read_dict[identifier][0])
    elif score1 < score2:
        print "identifier: " + identifier + "\t\test. best position: " + \
        str(best_position2) + "\t\tactual position: " + str(read_dict[identifier][0])
    elif score1 == score2 and score1 == -1:
        print "identifier: " + identifier + "\t\test. best position: " + \
        "NO MATCH" + "\t\tactual position: " + str(read_dict[identifier][0])                
    else:
        print "identifier: " + identifier + "\t\test. best position: " + \
        "MESSED UP" + "\t\tactual position: " + str(read_dict[identifier][0])
    """


def best_match_position(bw, bwr, s, diff, sa):
    sa_index_list = inexact_search(bw, bwr, s, diff)
    if len(sa_index_list) != 0:
        best_index, score = sa_index_list[0]
        return sa[best_index]+1, score
    else:
        return -1,-1
    

def main():

    align_reads(sys.argv[1], sys.argv[2])
    
main()

"""
def parse_genome():

    # Remove all newline characters
    if sys.argv[2] == "-c":

        file_name = sys.argv[1]
        new_file_name = file_name[:-4] + '_new.fna'
        read_file = open(file_name, 'r')
        write_file = open(new_file_name, 'w')

        for line in read_file:
            write_file.write(line.rstrip())

        read_file.close()
        write_file.close()

    # Remove metadata, only get reads
    elif sys.argv[2] == "-b":

        file_name = sys.argv[1]
        new_file_name = file_name[:-6] + '_dna.fastq'
        read_file = open(file_name, 'r')
        write_file = open(new_file_name, 'w')

        line_num = 0
        for line in read_file:
            line_num += 1
            if line_num % 4 == 2:
                write_file.write(line)

        read_file.close()
        write_file.close()
"""