# Chris Frey
# 10/8/2020
# Developed in Python 3.5.3
# Needleman-Wunsch Algorithm

import sys
import os

def read_seq(f_name):
    """
    Reads in both base or protein input sequences from passed file name.
    Ensures that both sequences are in uppercase form.
    Returns a tuple with both sequences in list form.
    
    Required file format: <row sequence>
                          <empty line>
                          <column sequence>

    Args: str f_name - Name of sequence input file

    Return: tuple(list row_seq, list col_seq) - Tuple of input sequences in list form.
    """

    with open(f_name, 'r') as seq_f:
        row_seq = seq_f.readline().strip()
        seq_f.readline()
        col_seq = seq_f.readline().strip()

    if not row_seq.isupper():
        row_seq = row_seq.upper()

    if not col_seq.isupper():
        col_seq = col_seq.upper()

    return (list(row_seq), list(col_seq))

def read_substn_matrix(f_name):
    """
    Reads in base or protein substitution matrix from passed file name.
    Ensures that base/protein references are in uppercase form.
    Returns a 2D list representing the substitution matrix.
    Row #1 and column #1 are considered to be equivalent, thus
    column #1 is discarded during the matrix read.
    
    Required file format:  Bases/Proteins
                          B val val val ...
                          a val val val ... 
                          s val val val ...
                          e val val val ...
                          / val val val ...
                          P val val val ...
                          r val val val ...
                          o val val val ...

    Args: str f_name - Name of substitution matrix input file

    Return: list(list) substn_mtrx - Substitution matrix in 2D list form.
    """

    substn_mtrx = []

    with open(f_name, 'r') as substn_f:
        frst_line = substn_f.readline().strip()

        if not frst_line.isupper():
            frst_line = frst_line.upper()
        
        substn_mtrx.append(frst_line.split())

        for line in substn_f:
            line = line.strip()

            if not line.isupper():
                line = line.upper()

            if line:
                substn_mtrx.append(line.split()[1:])

    return substn_mtrx

def build_substn_dict(substn_mtrx):
    """
    Converts the passed base or protein substitution matrix to 
    dictionary form and determines if the sequence is DNA or protein.
    Keys are choose-two permutations of the base or protein concatenations.
    Values are the substitution score for their associated key determined by the matrix.
    Returns a tuple containing the substitution dictionary and the sequence type.

    Args: list(list) substn_mtrx - Substitution matrix.

    Return: tuple(dict substn_dct, str seq_type) - Tuple of substitution dictionary and sequence type string.
    """

    substn_dct = {}

    for i in range(len(substn_mtrx[0])):
        for j in range(i+1, len(substn_mtrx)):
            substn_dct[substn_mtrx[0][i] + substn_mtrx[0][j-1]] = int(substn_mtrx[j][i])
            substn_dct[substn_mtrx[0][j-1] + substn_mtrx[0][i]] = int(substn_mtrx[j][i])

    seq_type = 'Protein' if len(substn_dct) == 400 else 'DNA'

    return (substn_dct, seq_type)

def build_F(row_seq, col_seq, g):
    """
    Builds the 'F' matrix and computes/fills in the first row 
    and column of alignment scoring values.
    Returns a 2D list representing the 'F' matrix.

    Args: list row_seq - Row sequence in list form.
          list col_seq - Column sequence in list form.
          int g - Gap penalty value.

    Return: list(list) F - 2D list representing the 'F' matrix for the NW algorithm.
    """
    
    col_seq.insert(0, ' ')
    col_seq.insert(0, ' ')
    F = [col_seq]
    row_seq.insert(0, ' ')
    
    for base in row_seq:
        F.append([0 for i in range(len(col_seq))])
        F[-1][0] = base

    for i in range(2, len(F[1])):
        F[1][i] = g * (i-1)

    for i in range(2, len(F)):
        F[i][1] = g * (i-1)

    return F

def _build_P(F):
    """
    Creates matrix 'P' for tracking path direction during NW algorithm 
    and computes/fills the first row and column of path direction values.
    Returns the matrix 'P' of size len(F)-1 on both dimensions.

    Args: list(list) F - F matrix in 2D list form.

    Return: list(list) P - P matrix in 2D list form.
    """

    P = [['l' for i in range(len(F[0])-1)]]
    P[0][0] = ''

    for i in range(len(F)-2):
        P.append(['' for i in range(len(F[0])-1)])
        P[-1][0] = 'u'
    
    return P

def _nw_fill(F, P, substn_dct, g):
    """
    In-place value computation and population of matrices F and P
    per the Needleman-Wunsch highroad alignment algorithm.

    Args: list(list) F - F matrix in 2D list form.
          list(list) P - P matrix in 2D list form.
          dict substn_dct - Substitution dictionary.
          g - Gap penalty value.
    """

    for i in range(2, len(F)):
        for j in range(2, len(F[0])):
            up = F[i-1][j] + g
            diag = F[i-1][j-1] + substn_dct[F[0][j] + F[i][0]]
            left = F[i][j-1] + g
            max_val = max(up, max(diag, left))

            if up == max_val:
                P[i-1][j-1] = 'u'
            elif diag == max_val:
                P[i-1][j-1] = 'd'
            else:
                P[i-1][j-1] = 'l'

            F[i][j] = max_val

def _nw_backtrace(P, F):
    """
    Performs backtrace of optimal alignment based on the bottom-
    right to top-left path in matrix P.
    Returns a tuple with alignment score, both sequences in alignment form,
    and a sequence denoting matching symbols.

    Args: list(list) - Matrix P in 2D list form.

    Return: tuple(int score, str x, str y, str m) - Tuple of alignment score, aligned sequences,
                                                    and match sequence.
    """
    i = len(P) - 1
    j = len(P[0]) - 1
    x, y, m = [], [], []
    score = F[i+1][j+1]

    while P[i][j] != '':
        # Swap first conditional block with trailing else
        # block (P[i][j] == 'l') for lowroad alignment.
        if P[i][j] == 'u':
            x.insert(0, '-')
            y.insert(0, F[i+1][0])
            m.insert(0, ' ')
            i -= 1
            
        elif P[i][j] == 'd':
            x.insert(0, F[0][j+1])
            y.insert(0, F[i+1][0])
            m.insert(0, '|') if x[0] == y[0] else m.insert(0, ' ')
            i -= 1
            j -= 1

        else:
            x.insert(0, F[0][j+1])
            y.insert(0, '-')
            m.insert(0, ' ')
            j -= 1
    
    return (score, ''.join(x), ''.join(y), ''.join(m))

def needleman_wunsch(F, substn_dct, g):
    """
    Wrapper function that calls/performs the major steps
    of the Needleman-Wunsch highroad alignment algorithm.
    Returns a tuple with alignment score, both sequences in alignment form,
    and a sequence denoting matching symbols.

    Args: list(list) F - F matrix in 2D list form.
          dict substn_dct - Substitution dictionary.
          g - Gap penalty value.

    Return: tuple(int score, str x, str y, str m) - Tuple of alignment score, aligned sequences,
                                                    and match sequence.
    """

    P = _build_P(F)
    _nw_fill(F, P, substn_dct, g)
    return _nw_backtrace(P, F)

def write_output(score, x, y, m, g, seq_type):
    """
    Writes optimal alignment sequences, symbol matches, and score to output text file.
    
    Out file format: Global alignment
                     DNA or Protein Sequence
                     g = <user defined value>
                     <column sequence>
                     <matching symbol sequence>
                     <row sequence>
                     Alignment score = <NM highroad alignment score>

    Args: int score - NM highroad alignment score.
          str x - Aligned column sequence.
          str y - Aligned row sequence.
          str m - Matching symbol sequence.
          int g - Gap penalty value.
          str seq_type - Sequence type (DNA or Protein).
    """

    with open('../Results/Human-mouse-alightment.txt', 'w+') as out_f:
            out_f.write('Global alignment\n')
            out_f.write(seq_type + ' Sequence\n')
            out_f.write('Gap penalty = ' + str(g) + '\n')
            out_f.write(x + '\n')
            out_f.write(m + '\n')
            out_f.write(y + '\n')
            out_f.write('Alignment score = ' + str(score))

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print('Error: Unknown program argument form.')
        print('Required Form: <program_file.py> <sequence_file.txt> <sub_matrix.txt> <gap penalty value>')
        sys.exit(0)

    elif not os.path.isfile(sys.argv[1]):
        print('Error: Unable to find sequence text file. Please check program arguments.')
        sys.exit(0)

    elif not os.path.isfile(sys.argv[2]):
        print('Error: Unable to find substitution matrix text file. Please check program arguments.')
        sys.exit(0)

    g = 0

    try:
        g = int(sys.argv[3])

    except:
        print('Error: Invalid gap penalty value. Please check program arguments.')
        sys.exit(0)
    
    seq_fn = sys.argv[1]
    substn_mtrx_fn = sys.argv[2]

    row_seq, col_seq = read_seq(seq_fn)
    F = build_F(row_seq, col_seq, g)
    substn_mtrx = read_substn_matrix(substn_mtrx_fn)
    substn_dct, seq_type = build_substn_dict(substn_mtrx)
    score, x, y, m = needleman_wunsch(F, substn_dct, g)
    write_output(score, x, y, m, g, seq_type)
