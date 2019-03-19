import sys

## Local alignment implementation but,
## only finds aligned sequence with the highest score
## # TODO: Find other max values to use for local alignments ????
##         Threshold for what a max value can be????
##         FROM wikipedia:
##         To obtain the second best local alignment, apply the traceback process starting at
##         the second highest score outside the trace of the best alignment.
## # TODO: Global alignment?????
##         FROM wikipedia: https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
##         Probably easy, just need to create another matrix then traceback from,
##         absolute bottom right.
## # TODO: Add functions for adding statistics to final output;
##         how many similar 'chars', how many gaps, how many mismatches
##         total score

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
    # Some sequence info
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

    ## Save this for now for debugging purposes
    #prints matrix all clean like
    #print('        ', ''.join(['{:5}'.format(item) for item in seq2]))
    #print('\n'.join([''.join(['{:5}'.format(item) for item in row]) for row in matrix]))

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
    ##  if stepper   => diagonal: match/mismatch
    ##  elif stepper => up: gap in seq1
    ##  elif stepper => left: gap in seq2

    align1 = []
    align2 = []
    (x,y) = maxpos
    stepper = next_step(matrix, x, y)

    # if our max step is diagonally
    #   from seq1,seq2 append both 'chars' to align1, and align2
    # if our max step is up
    #   append a gap to align1
    #   from seq2 append 'char' to align2
    # if our max step is left
    #   from seq1 append 'char' to align1
    #   append a gap to align2
    while(stepper != 0):
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
        stepper = next_step(matrix, x, y)

    # Appends last 'chars' to our align1,align2 after
    # while loop is dropped out of by the return of 0 from next_step()
    align1.append(seq1[x-1])
    align2.append(seq2[y-1])

    # Initially reversed because we are working backward on the string
    # Use array manipulation to reverse string
    align1 = align1[::-1]
    align2 = align2[::-1]

    # Creates a 'formatter' for printing out pretty comparison.
    formatter = []
    for i in range(len(align1)):
        if align1[i] == align2[i]:
            formatter.append('|')
        elif (align1[i] == '-') or (align2[i] == '-'):
            formatter.append(' ')
        elif align1[i] != align2[i]:
            formatter.append(':')

    ###### DIRTY(DON'T LOOK) ######
    ## Move this to another function write_out() ???? Probably
    ## Also do this all in one go using maxLen and by20
    maxLent = len(align1)
    by20t = 0
    while by20t < maxLent/20:
        print(by20t*20 + 1, ':', ((by20t+1)*20)-1)
        print('   ', end='')
        for s in align1[20*by20t:by20t*20 +19]:
            print(s, end='')
        print('\n')
        print('   ', end='')
        for t in formatter[20*by20t:by20t*20 +19]:
            print(t, end='')
        print('\n')
        print('   ', end='')
        for u in align2[20*by20t:by20t*20 +19]:
            print(u, end='')
        print('\n')

        by20t += 1

    maxLen = len(align1)
    by20 = 0
    ###### MORE DIRT FOR THE WRITE FILE######
    with open("segment_aligned.txt", 'w+') as f:
        while by20 < maxLen/20:
            f.write(str(by20*20 +1))
            f.write(':')
            f.write(str(((by20+1)*20) -1))
            f.write('\n')
            for s in align1[20*by20:by20*20 +19]:
                f.write(s)
            f.write('\n')
            for t in formatter[20*by20:by20*20 +19]:
                f.write(t)
            f.write('\n')
            for u in align2[20*by20:by20*20 +19]:
                f.write(u)
            f.write('\n\n')
            by20 += 1
    return 0

def next_step(matrix, x, y):
    diag =  matrix[x-1][y-1]
    up   =  matrix[x][y-1]
    left =  matrix[x-1][y]
    #print(diag, ' ', up, ' ', left)
    maxstep = max(diag,up,left)
    # Quick dropout function. If diag or up or left is zero
    # then the function ends traceback.
    if (diag and up and left) == 0:
        print('Finished with traceback')
        return 0

    # where diag is matrix[x-1][y-1]
    #       up   is matrix[x][y-1]
    #       left is matrix[x-1][y]
    #       maxstep is whatever is the maximum between the 3 (diag,up,left)
    # if diag == maxstep -> appends both 'chars'
    #                       from seq1 and seq2 to align1 and align2
    # if up   == maxstep -> appends a gap for align1
    #                       and appends 'char' from seq2 to align2
    # if left == maxstep -> appends 'char' from seq1 to align1
    #                       and appends a gap for align2
    if diag == maxstep:
        return 1
    if up == maxstep:
        return 2
    if left == maxstep:
        return 3
    assert False

def score_matrix(matrix, x, y):
    # Scores the matrix. Starts at top left of matrix and works its way
    # right until end of column. Then continues to the next row.
    # Returns the max of the diagonal value, up_adj, and left_adj
    # (given match, mismatch, and gap penalty)

    # if seq1[a] == seq2[b]
    #   then: match
    #   else: mismatch
    diag     = matrix[x-1][y-1] + (match if seq1[x-1] == seq2[y-1] else mismatch)

    ## These are flipped on accident (i think????)
    up_adj   = matrix[x-1][y] + gap
    left_adj = matrix[x][y-1] + gap
    return max(0, diag, up_adj, left_adj)

def read_file():
    # Reading stuff from file. Has to be one long string
    # in text file.
    ## Use delimeter for multi-line reading????
    f = open("segment.txt", "r")
    x = f.read().split('\n')
    Seq1 = x[0]
    Seq2 = x[1]

    return Seq1, Seq2

if __name__ == "__main__":

    seq1 = 'tgttacgg'
    seq2 = 'ggttgacta'
    #seq1, seq2 = read_file()

    matrix, maxpos = align(seq1, seq2)
    traceback(matrix, maxpos)

    print('Match if: \'|\'')
    print('Mismatch if \':\'')
    print('Gap if space')
