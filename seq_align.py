import sys
# For my comments, usually one # (#) means comment is to explain code,
# and two # (##) non-permanent stuff

##  SHOULD PROBABLY ADD TESTS TO CONFIRM THIS EVEN works

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
##         how many similar bases, how many gaps, how many mismatches
##         total score
## # TODO: Rename align(), traceback() so that global alignment can be implemented,
##         align()-> localAlign(), traceback()-> localTraceback(), score_matrix() -> localMatrix()
##                -> globalAlign(),           -> globalTraceback,                 -> globalMatrix()
##         Make scores a local variable. (Use same scores for local and global???)
##         Changeable scores???

#match      = +m
localmatch       = 3
#mismatch   = -s
localmismatch    = -3
#gap        = -d
localgap         = -2
#scoring  F = (# matches)*m - (# mismatches)*s -(#gaps)*d

globalmatch = 1
globalmismatch = -1
globalgap = -1

##minimal edit distance: given two strings x, y, find minimum # of edits
##(insertions,deletions,mutations) to transform one string to the other

##ok to have unlimited # of gaps in beginning & end?
##don't penalize gaps on either end

##insertions and/or deletions are called indels

def local_align(seq1, seq2):
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

    # NOTE: 1: Determine gap penalty scheme.
    #          * s(a,b) - similarity score of the elements that
    #            constituted the two sequences
    #          * W_k - the penalty of a gap that has length k
    # GLOBAL VARIABLES AT THE TOP
    # match       = 3
    # mismatch    = -3
    # gap         = -2

    # NOTE: 2: & 3: creates a scoring matrix and fills according to
    #               match, mismatch, and gap scoring. Checks the max
    #               from the diagonal, left_adj, and up_adj. If none
    #               are above 0 fill current spot in matrix with 0.
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
    ##prints matrix all clean like
    print('        ', ''.join(['{:5}'.format(item) for item in seq2]))
    print('\n'.join([''.join(['{:5}'.format(item) for item in row]) for row in matrix]))

    return matrix, max_pos, max_score

def local_traceback(matrix, maxpos):
    # NOTE: 4: Traceback in the matrix to find best alignment.
    #          Start at element with highest score, traceback based on the
    #          source of each score recursively, until 0 is encountered. The
    #          segments that have the highest similarity score based on the
    #          given scoring system is generated in this process. To obtain
    #          the second best local alignment, apply the traceback process
    #          starting at the second highest score outside the trace of the
    #          best alignment.
    #
    # Stepper looks for max(up,left,diagonal)
    #
    # while(0 is not found)
    #  if stepper   => diagonal: match/mismatch
    #  elif stepper => up: gap in seq1
    #  elif stepper => left: gap in seq2

    align1 = []
    align2 = []
    (x,y) = maxpos
    stepper = next_step(matrix, x, y)

    # if our max step is diagonally
    #   from seq1,seq2 append both bases to align1, and align2
    # if our max step is up
    #   append a gap to align1
    #   from seq2 append base to align2
    # if our max step is left
    #   from seq1 append base to align1
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

    # Appends last bases to our align1,align2 after
    # while loop is dropped out of by the return of 0 from next_step()
    align1.append(seq1[x-1])
    align2.append(seq2[y-1])

    # Initially reversed because we are working backward on the string
    # Use array manipulation to reverse string
    align1 = align1[::-1]
    align2 = align2[::-1]

    return align1, align2

# Takes the next step in the traceback
def next_step(matrix, x, y):
    diag =  matrix[x-1][y-1]
    up   =  matrix[x][y-1]
    left =  matrix[x-1][y]

    ## WHAT IF TWO ARE THE SAME???
    ## If diag is the same with either(up, left) always go diag
    ## What if up is same as left??? With max() I think hierarchy is diag->up->left
    ## Does it matter???
    maxstep = max(diag,up,left)
    # Quick reject function. If diag or up or left is zero
    # then the function ends traceback.
    if (diag and up and left) == 0:
        print('Finished with traceback')
        return 0

    # where diag is matrix[x-1][y-1]
    #       up   is matrix[x][y-1]
    #       left is matrix[x-1][y]
    #       maxstep is whatever is the maximum between the 3 (diag,up,left)
    # if diag == maxstep -> appends both bases
    #                       from seq1 and seq2 to align1 and align2
    # if up   == maxstep -> appends a gap for align1
    #                       and appends base from seq2 to align2
    # if left == maxstep -> appends base from seq1 to align1
    #                       and appends a gap for align2
    if diag == maxstep:
        return 1
    if up == maxstep:
        return 2
    if left == maxstep:
        return 3
    assert False

# Scores matrix for local alignment
def score_matrix(matrix, x, y):
    # Scores the matrix. Starts at top left of matrix and works its way
    # right until end of column. Then continues to the next row.
    # Returns the max of the diagonal value, up_adj, and left_adj
    # (given match, mismatch, and gap penalty)

    # if seq1[a] == seq2[b]
    #   then: match
    #   else: mismatch
    diag     = matrix[x-1][y-1] + (localmatch if seq1[x-1] == seq2[y-1] else localmismatch)

    ## These are flipped on accident (i think????)
    up_adj   = matrix[x-1][y] + localgap
    left_adj = matrix[x][y-1] + localgap
    return max(0, diag, up_adj, left_adj)

# Reading stuff from file. Has to be one long string in text file.
def read_file():
    ## Use delimeter for multi-line reading????
    f = open("segment.txt", "r")
    x = f.read().split('\n')
    Seq1 = x[0]
    Seq2 = x[1]

    return Seq1, Seq2

# Write aligned sequence to console and to text file (segment_aligned.txt)
def outputs(align1, align2):
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
    ## There's gotta be a better way......
    ## Move this to another function write_out() ???? Probably
    ## Also do this all in one go using maxLen and by20. maybe not.
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

def scores(align1, align2, max_score):
    ##scores; matches, mismatches, indels
    totalScore = 0
    same = 0

    for i in range(0,len(align1)):
        if align1[i] == align2[i]:
            same = same + 1
            totalScore = totalScore + localmatch
        elif (align1[i] != align2[i]) and ((align1[i] and align2[i]) != '-'):
            totalScore = totalScore + localmismatch
        elif align1[i] or align2[i] == '-':
            totalScore = totalScore + localgap
    identity = float(same) / len(align1) * 100
    print('--Smith-Waterman Results--')
    print('Total Score:', totalScore)
    print('max_score:', max_score)
    print('Identical bases: ', same)
    print('Percent Identity: %.2f' % identity, '%')
    print('Aligned Sequence Length: ', len(align1))
    print('\n')

def globalAlign(seq1, seq2):
    seq1 = "-" + seq1
    seq2 = "-" + seq2
    seq1len = len(seq1)
    seq2len = len(seq2)
    matrix = []
    # Needleman-Wunsch algorithm for global sequence alignment
    # https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
    # 1: Construct grid using seq1 and seq2
    # 1.5: choose scoring system,
    #      letters may either match, mismatch, or be matched to a gap.
    #      Let's use Match = +1, Mismatch or Indel: -1
    # 2: Filling in the table. First row and first column are scored by using
    #    mismatch penalty increasing(decreasing) by what the penalty is
    #    Each cell has 3 possible candidates.
    #      1: Diagonal top-left neighbor has (x-1,y-1). Then comparing the two
    #         bases. Add match if bases are same to top-left neighbor
    #      2: Top-neighbor has a score (x, y-1). Moving from there represents an
    #         indel so add mismatch score to top-neighbor
    #      3: Left-neight has a score (x-1, y). Moving from there represents an
    #         indel so add mismatch score to left-neighbor
    #    Highest score is entered into the cell.
    #    IF(highest candidate score is from two or all neighboring cells)
    #      THEN: All directions reaching the highest candidate score must be noted
    #            as possible origin cells in the finished matrix.
    # 3: Traceback, tracing arrows back to origin,
    #    Make a path from the cell on the bottom right bacl to the cell on the top left
    #    by following the direction of the arrows,
    #    Sequence is constructed by these rules:
    #       1: A diagonal arrow represents a match or mismatch so the letter of the
    #          column and the letter of the row of the origin cell will align
    #       2: a horizontal or vertical arrow represents an indel. Horizonal will
    #          align a gam ('-') to the letter of the row, vertical arrows will align
    #          a gap to the letter of the column.
    #       3: If there are multiple arrows to choose from, they represent a branching of
    #          the alignments. If two or more branches all belong to paths from the bottom
    #          left to the top right cell, they are equally viable alignments.
    for i in range(seq1len):
         matrix.append([0]*seq2len)

    # Fills in first row and first column with necessary score
    matrix[0][0] = 0
    for i in range(1,seq1len):
        matrix[0][i] = matrix[0][i-1] + globalgap
    for j in range(1,seq2len):
        matrix[j][0] = matrix[j-1][0] + globalgap

    # Scores the rest of the table
    for i in range(1, seq1len):
        for j in range(1, seq2len):
            diag        = matrix[i-1][j-1] + globalmatch
            left_adj    = matrix[i-1][j] + globalgap
            top_adj     = matrix[i][j-1] + globalgap
            matrix[i][j] = max(diag, left_adj, top_adj)
    return matrix

def g_next_step(matrix, x, y):
    current = matrix[x][y]
    diag = matrix[x-1][y-1]
    up   = matrix[x][y-1]
    left = matrix[x-1][y]

    if(x==0 or y==0):
        return 0
    print('current', current)
    print('diag', diag + globalmatch)
    print('up', up -globalgap)
    print('left',left-globalgap)

    if current == (diag + globalmatch):
        return 1
    elif current == (up - globalgap):
        return 2
    elif current == (left - globalgap):
        return 3

    assert False

def globalTraceback(matrix, seq1, seq2):
    align1 = []
    align2 = []
    seq1len = len(seq1)
    seq2len = len(seq2)
    x = seq1len - 1
    y = seq2len - 1

    step = g_next_step(matrix, x, y)
    while(step!= 0):
        if step == 1:
            align1.append(seq1[x-1])
            align2.append(seq2[y-1])
            x-=1
            y-=1
        elif step == 2:
            print('-')
            align1.append('-')
            align2.append(seq2[y-1])
            y-=1
        elif step == 3:
            print('-')
            align1.append(seq1[x-1])
            align2.append('-')
            x-=1
        step = g_next_step(matrix, x, y)


    align1 = align1[::-1]
    align2 = align2[::-1]

    align1.append(seq1[x-1])
    align2.append(seq2[y-1])

    return align1, align2

if __name__ == "__main__":

    ##seq1 = 'tgttacgg'
    ##seq2 = 'ggttgacta'
    seq1, seq2 = read_file()

    # matrix, maxpos, max_score = local_align(seq1, seq2)
    # align1, align2 = local_traceback(matrix, maxpos)
    # outputs(align1, align2)
    # scores(align1,align2, max_score)

    seq1 = 'GCATGCU'
    seq2 = 'GATTACA'
    matrix = globalAlign(seq1, seq2)

    print('        ', ''.join(['{:5}'.format(item) for item in seq2]))
    print('\n'.join([''.join(['{:5}'.format(item) for item in row]) for row in matrix]))

    align1, align2 = globalTraceback(matrix,seq1,seq2)

    print(align1, '\n', align2)
