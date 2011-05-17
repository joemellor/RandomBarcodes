__author__ = 'joemellor'

import random
import itertools
import Bio
from Bio import pairwise2

def hamming_distance(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))
def dameraulevenshtein_distance(seq1, seq2):
    # codesnippet:D0DE4716-B6E6-4161-9219-2903BF8F547F
    # Conceptually, this is based on a len(seq1) + 1 * len(seq2) + 1 matrix.
    # However, only the current and two previous rows are needed at once,
    # so we only store those.
    oneago = None
    thisrow = range(1, len(seq2) + 1) + [0]
    for x in xrange(len(seq1)):
        # Python lists wrap around for negative indices, so put the
        # leftmost column at the *end* of the list. This matches with
        # the zero-indexed strings and saves extra calculation.
        twoago, oneago, thisrow = oneago, thisrow, [0] * len(seq2) + [x + 1]
        for y in xrange(len(seq2)):
            delcost = oneago[y] + 1
            addcost = thisrow[y - 1] + 1
            subcost = oneago[y - 1] + (seq1[x] != seq2[y])
            thisrow[y] = min(delcost, addcost, subcost)
            # This block deals with transpositions
    return thisrow[len(seq2) - 1]

barcodelist1 = list()
f = open('300barcodes.txt', 'r')
for line in f:
    barcodelist1.append(line.rstrip("\n"))
f.close()
barcodelist2 = barcodelist1
f = open('temp2.txt', 'w')
print len(barcodelist1)
f.write("\t")
for s1 in barcodelist1:
    f.write(s1)
    f.write("\t")
f.write("\n")
for s1 in barcodelist1:
    f.write(s1)
    f.write("\t")
    for s2 in barcodelist2:
        align = pairwise2.align.globalms(s1,s2,1,-1,-1,-1)
        seq1, seq2, score, start, end = align[0]
        f.write(str(score))
        f.write("\t")
    f.write("\n")
f.close()

