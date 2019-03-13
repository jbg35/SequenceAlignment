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

#insertions and/or deletions are called indels



def align(seq1, seq2):
    print('seq1: ', seq1,'\n' 'seq2: ',  seq2)
    seq1len = len(seq1)
    seq2len = len(seq2)
    print('len: ', seq1len, seq2len)

    #initialize matrix for scoring
    ## FIXME: make size for only what is needed????
    matrix = [[0 for x in range(max(seq1len+1,seq2len+1))] for y in range(max(seq1len+1,seq2len+1))]
    for i in range (1, seq1len+1):
        matrix[i][0] = i * gap
        matrix[0][i] = i * gap

    for i in range(seq1len+1):
        for j in range(seq2len+1):
            if seq1[i-1] == seq2[j-1]:
                matrix[i][j] = matrix[i-1][j-1]
            else:
                matrix[i][j] = min({matrix[i-1][j-1]+mismatch,
                                    matrix[i-1][j]+ gap      ,
                                    matrix[i][j-1]+ gap       })

    #prints matrix all clean like
    print('   ', ''.join(['{:5}'.format(item) for item in seq1]))
    print('\n'.join([''.join(['{:5}'.format(item) for item in row]) for row in matrix]))
    print('\n'.join('{:}'.format(test) for test in seq1 ) )

if __name__ == "__main__":

    seq1 = 'atgagtcttc'
    seq2 = 'atcgacgtca'

    align(seq1, seq2)
