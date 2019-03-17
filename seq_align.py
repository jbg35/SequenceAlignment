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
    #Printing string info
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
    #    The max of the diagonals, left adjacent and right adjacent, or 0
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
    print(matrix[7][6])
    #prints matrix all clean like
    print('        ', ''.join(['{:5}'.format(item) for item in seq2]))
    ## TODO: left letters maybe???
    print('\n'.join([''.join(['{:5}'.format(item) for item in row]) for row in matrix]))

    return matrix, max_pos

def traceback(matrix, maxpos):
    ## NOTE: 4: Traceback in the matrix to find best alignment.
    ##          Start at element with highest score, traceback based on the
    ##          source of each score recursively, until 0 is encountered. The
    ##          segments that have the highest similarity score based on the
    ##          given scoring system is generated in this process. To obtain
    ##          the second best local alignment, apply the traceback process
    ##          starting at the second highest score outside the trace of the
    ##          best alignment.
    ##
    ## Stepper looks for max(up,left,diagonal)
    ##
    ## while(0 is not found)
    ##  if stepper => diagonal: match/mismatch
    ##  elif stepper => up: gap in seq1
    ##  elif stepper => left: gap in seq2
    align1 = []
    align2 = []
    (x,y) = maxpos
    stepper = next_step(matrix, x, y)
    while(stepper != 0):
        # Keep char in string if we traverse diagonaly
        print('step:', stepper)
        if stepper == 1:
            align1.append(seq1[x-1])
            align2.append(seq2[y-1])
            x -= 1
            y -= 1
        elif stepper == 2:
            align1.append('-')
            align2.append(seq2[y-1])
            y -= 1
        elif stepper == 3:
            align1.append(seq1[x-1])
            align2.append('-')
            x -= 1
        # print(align1)
        # print()
        # print(align2)
        stepper = next_step(matrix, x, y)

        print('step after:', stepper)


    align1.append(seq1[x-1])
    align2.append(seq2[y-1])
    # Initially reversed because we are working backward on the string
    print(align1)
    print(align2)
    print('\n')

    align1 = align1[::-1]
    align2 = align2[::-1]

    formatter = []
    for i in range(len(align1)):
        if align1[i] == align2[i]:
            formatter.append('|')
        elif (align1[i] or align2[i]) == '-':
            formatter.append(' ')
        elif align1[i] != align2[i]:
            formatter.append(':')

    print(' '.join(align1))
    print(' '.join(formatter))
    print(' '.join(align2))

    return 0

def next_step(matrix, x, y):
    diag =  matrix[x-1][y-1]
    up   =  matrix[x][y-1]
    left =  matrix[x-1][y]
    print(diag, ' ', up, ' ', left)
    maxstep = max(diag,up,left)
    if (diag and up and left) == 0:
        print('END')
        return 0

    if diag == maxstep:
        return 1 
    if up == maxstep:
        return 2 
    if left == maxstep:
        return 3 
    assert False

def score_matrix(matrix, x, y):
    # Scores the matrix. Starts at top left of matrixand works its way
    # right until end of column. Then continues to the next row.
    # Returns the max of the diagonal value, up_adj, and left_adj
    # (given match, mismatch, and gap penalty)

    # if seq1[a] == seq2[b]
    #   then: match
    #   else: mismatch
    diag     = matrix[x-1][y-1] + (match if seq1[x-1] == seq2[y-1] else mismatch)
    up_adj   = matrix[x-1][y] + gap
    left_adj = matrix[x][y-1] + gap
    return max(0, diag, up_adj, left_adj)

if __name__ == "__main__":

    seq1 = 'tgttacgg'
    seq2 = 'ggttgacta'

    matrix, maxpos = align(seq1, seq2)
    traceback(matrix, maxpos)
