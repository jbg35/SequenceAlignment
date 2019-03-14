import sys
#match      = +m
match       = 3
#mismatch   = -s
mismatch    = -3
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


    # Smith-Waterman algorithm for local sequence alignment
    # https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
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
    # 4: Traceback in the matrix. starting at the element with the highest
    #    score, traceback based on the source of each score reursively,
    #    until 0 is encountered

    ## NOTE: 1: Determine gap penalty scheme.
    ##          * s(a,b) - similarity score of the elements that
    ##            constituted the two sequences
    ##          * W_k - the penalty of a gap that has length k
    # GLOBAL VARIABLES AT THE TOP
    # match       = 3
    # mismatch    = -3
    # gap         = -2

    ## NOTE: 2: & 3: creates a scoring matrix and fills according to
    ##               match, mismatch, and gap scoring. Checks the max
    ##               from the diagonal, left_adj, and up_adj. If none
    ##               are above 0 fill current spot in matrix with 0.
    matrix = [[0 for x in range(seq2len+1)] for y in range(seq1len+1)]
    max_score = 0
    max_pos = None
    for i in range(1, seq1len+1):
        for j in range(1,seq2len+1):
            score  = score_matrix(matrix,i,j)
            if score > max_score:
                max_score = score
                max_pos = (i,j)
            matrix[i][j] = score

    ## NOTE: 4: Traceback in the matrix to find best alignment.
    ##          Start at element with highest score, traceback based on the
    ##          source of each score recursively, until 0 is encountered. The
    ##          segments that have the highest similarity score based on the
    ##          given scoring system is generated in this process. To obtain
    ##          the second best local alignment, apply the traceback process
    ##          starting at the second highest score outside the trace of the
    ##          best alignment.
    ## TODO: TRACEBACK


    #prints matrix all clean like
    print('        ', ''.join(['{:5}'.format(item) for item in seq2]))
    ## TODO: Left side of matrix letters.
    print('\n'.join([''.join(['{:5}'.format(item) for item in row]) for row in matrix]))

    return matrix, max_pos

def score_matrix(matrix, x, y):
    similarity = match if seq1[x-1] == seq2[y-1] else mismatch
    diag = matrix[x-1][y-1] + similarity
    up_adj = matrix[x-1][y] + gap
    left_adj = matrix[x][y-1] + gap
    return max(0, diag, up_adj, left_adj)

if __name__ == "__main__":

    seq1 = 'ggttgacta'
    seq2 = 'tgttacgg'

    align(seq1, seq2)
