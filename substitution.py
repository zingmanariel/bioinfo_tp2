"""Substitution matrices for the protein alignments (3.2 a).

Biopython ships the matrices already loaded, so there's no need to download
them:

    from Bio.Align import substitution_matrices
    substitution_matrices.load()          # available names

Note: the TP suggests "Blossum60", which does not exist among the standard
ones. The BLOSUM matrices available are 45, 50, 62, 80 and 90; we default to
BLOSUM62, which is what BLASTP uses and therefore what makes our results
comparable to BLAST's. The number is the identity percentage at which the
sequences used to build the matrix were clustered: lower (BLOSUM45) for
distant proteins, higher (BLOSUM80) for close ones.
"""

from Bio.Align import substitution_matrices

DEFAULT_MATRIX = 'BLOSUM62'

# Residue used to score symbols the matrix doesn't know.
UNKNOWN = 'X'


def available():
    """Names of the matrices Biopython ships."""
    return sorted(substitution_matrices.load())


def load_matrix(name=DEFAULT_MATRIX):
    """Loads the substitution matrix by name."""
    return substitution_matrices.load(name)


def score(matrix, a, b):
    """Score of substituting amino acid a for b.

    Symbols the matrix doesn't know are scored as 'X' (unknown residue),
    which is what real aligners do.
    """
    alphabet = matrix.alphabet
    a = a if a in alphabet else UNKNOWN
    b = b if b in alphabet else UNKNOWN
    return float(matrix[a, b])


def describe(matrix, name=DEFAULT_MATRIX):
    """Readable printout of the matrix, to look it over before using it."""
    print(f'Matrix {name}: alphabet of {len(matrix.alphabet)} symbols')
    print(matrix)

    diagonal = [matrix[a, a] for a in matrix.alphabet]
    off_diagonal = [matrix[a, b] for a in matrix.alphabet
                    for b in matrix.alphabet if a != b]
    print(f'\ndiagonal (same residue): min {min(diagonal):.0f}, '
          f'max {max(diagonal):.0f}, average {sum(diagonal) / len(diagonal):.2f}')
    print(f'off diagonal: min {min(off_diagonal):.0f}, '
          f'max {max(off_diagonal):.0f}, '
          f'average {sum(off_diagonal) / len(off_diagonal):.2f}')
    positives = sum(1 for value in off_diagonal if value > 0)
    print(f'distinct substitutions with a positive score (conservative): '
          f'{positives} of {len(off_diagonal)}')
