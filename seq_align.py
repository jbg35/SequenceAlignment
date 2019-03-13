import sys
#match      = +m
match       = 2
#mismatch   = -s
mismatch    = 0
#gap        = -d
gap         = -2
#scoring  F = (# matches)*m - (# mismatches)*s -(#gaps)*d

#minimal edit distance: given two strings x, y, find minimum # of edits
#(insertions,deletions,mutations) to transform one string to the other

#ok to have unlimited # of gaps in beginning & end?
#don't penalize gaps on either end
def align(seq1, seq2):
    print(seq1, seq2)
    seq1len = len(seq1)
    seq2len = len(seq2)
    print(seq1len, seq2len)

    #initialize matrix for scoring
    matrix = [[0 for x in range(seq1len)] for y in range(seq2len)]


    print(matrix[0][0])
    pass



if __name__ == "__main__":

    seq1 = 'atgagtcttc'
    seq2 = 'atcgacgtca'

    align(seq1, seq2)
