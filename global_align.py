import sys

globalmatch = 1
globalmismatch = -1
globalgap = -1

def globalAlign(seq1, seq2):
    # Adds gap for filling out matrix
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
    for j in range(1,seq2len):
        matrix[0][j] = matrix[0][j-1] + globalgap
    for i in range(1,seq1len):
        matrix[i][0] = matrix[i-1][0] + globalgap

    # Scores the rest of the table
    for i in range(1, seq1len):
        for j in range(1, seq2len):
            match     = matrix[i-1][j-1] + (globalmatch if seq1[i-1] == seq2[j-1] else globalmismatch)
            delete    = matrix[i-1][j] + globalgap
            insert     = matrix[i][j-1] + globalgap
            matrix[i][j] = max(match, delete, insert)
    return matrix

def g_next_step(matrix, x, y):
    current = matrix[x][y]
    diag = matrix[x-1][y-1]
    up   = matrix[x][y-1]
    left = matrix[x-1][y]

    if(x==0 or y==0):
        return 0
    # print('current', current)
    # print('diag', diag + globalmatch)
    # print('up', up -globalgap)
    # print('left',left-globalgap)

    if current == (diag + globalmatch):
        print((x,y), current, 'diag->', diag+ globalmatch)
        return 1
    elif current == (up + globalgap):
        print((x,y), current, 'up->', up+ globalgap)
        return 2
    elif current == (left + globalgap):
        print((x,y), current, 'left->', left+ globalgap)
        return 3
    else:
        print("what the",(x,y), current,diag + globalmatch, up + globalgap,left + globalgap)
    assert False

def globalTraceback(matrix, seq1, seq2):
    align1 = []
    align2 = []
    seq1len = len(seq1)
    seq2len = len(seq2)
    x = seq1len
    y = seq2len

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

    ## GLOBAL ALIGNMENT ##
    seq2 = 'GCATGCU'
    seq1 = 'GATTACA'
    matrix = globalAlign(seq1, seq2)
    print('        ', ''.join(['{:5}'.format(item) for item in seq2]))
    print('\n'.join([''.join(['{:5}'.format(item) for item in row]) for row in matrix]))
    align1, align2 = globalTraceback(matrix,seq1,seq2)
    print(align1, '\n', align2)
