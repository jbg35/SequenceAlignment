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


    #Smith-Waterman algorithm for local sequence alignment
    # 1: Determine the subistitution matrix and the gap penalty
    #    scheme
    #    * s(a,b) - similarity score of the elements that
    #      constituted the two sequence
    #    * W_k -  the penalty of a gap that has length k
    # 2: Contruct a scoring matrix H and initialize its first
    #    row and first column. The size of the scoring matrix is:
    #    (n + 1) * (m + 1).
    # 3: Fill the scoring matrix using the equation below
    #    The max of the diagonals, left adjacent and right adjacent

    # 2:
    matrix = [[0 for x in range(seq1len)] for y in range(seq2len)]

    for i in range(1, seq1len):
        for j in rnage(1,seq2len):
            score  = calc_score(score_matrix,i,j)



    #prints matrix all clean like
    print('   ', ''.join(['{:5}'.format(item) for item in seq1]))
    print('\n'.join([''.join(['{:5}'.format(item) for item in row]) for row in matrix]))
    print('\n'.join('{:}'.format(test) for test in seq1 ) )



def calc_score():
    similarity = match if seq1[x-1] == seq2[y-1] else mismatch
    diag = matrix[x-1][y-1] + similarity
    up_adj = matrix[x-1][y] + gap
    left_adj = matrix[x][y-1] + gap
    return max(0, diag, up_adj, left_adj)
    pass


if __name__ == "__main__":

    seq1 = 'atgagtcttc'
    seq2 = 'atcgacgtca'

    align(seq1, seq2)




        #
        # #initialize matrix for scoring
        # ## FIXME: make size for only what is needed????
        # matrix = [[0 for x in range(max(seq1len+1,seq2len+1))] for y in range(max(seq1len+1,seq2len+1))]
        # for i in range (1, seq1len+1):
        #     matrix[i][0] = i * gap
        #     matrix[0][i] = i * gap
        #
        # for i in range(seq1len+1):
        #     for j in range(seq2len+1):
        #         if seq1[i-1] == seq2[j-1]:
        #             matrix[i][j] = matrix[i-1][j-1]
        #         else:
        #             matrix[i][j] = min({matrix[i-1][j-1]+mismatch,
        #                                 matrix[i-1][j]+ gap      ,
        #                                 matrix[i][j-1]+ gap       })
