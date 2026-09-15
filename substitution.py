"""Matrices de sustitucion para los alineamientos de proteinas (3.2 a).

Biopython trae las matrices ya cargadas, asi que no hace falta bajarlas:

    from Bio.Align import substitution_matrices
    substitution_matrices.load()          # nombres disponibles

Nota: el TP sugiere "Blossum60", que no existe entre las estandar. Las BLOSUM
disponibles son 45, 50, 62, 80 y 90; usamos BLOSUM62 por defecto, que es la que
usa BLASTP y por lo tanto la que hace comparables nuestros resultados con los
de BLAST. El numero es el porcentaje de identidad al que se agruparon las
secuencias con las que se construyo la matriz: mas bajo (BLOSUM45) para
proteinas lejanas, mas alto (BLOSUM80) para proteinas cercanas.
"""

from Bio.Align import substitution_matrices

DEFAULT_MATRIX = 'BLOSUM62'

# Residuo con el que se puntuan los simbolos que la matriz no conoce.
UNKNOWN = 'X'


def available():
    """Nombres de las matrices que trae Biopython."""
    return sorted(substitution_matrices.load())


def load_matrix(name=DEFAULT_MATRIX):
    """Levanta la matriz de sustitucion por nombre."""
    return substitution_matrices.load(name)


def score(matrix, a, b):
    """Puntaje de sustituir el aminoacido a por el b.

    Los simbolos que la matriz no conoce se puntuan como 'X' (residuo
    desconocido), que es lo que hacen los alineadores.
    """
    alphabet = matrix.alphabet
    a = a if a in alphabet else UNKNOWN
    b = b if b in alphabet else UNKNOWN
    return float(matrix[a, b])


def describe(matrix, name=DEFAULT_MATRIX):
    """Impresion legible de la matriz, para analizarla antes de usarla."""
    print(f'Matriz {name}: alfabeto de {len(matrix.alphabet)} simbolos')
    print(matrix)

    diagonal = [matrix[a, a] for a in matrix.alphabet]
    off_diagonal = [matrix[a, b] for a in matrix.alphabet
                    for b in matrix.alphabet if a != b]
    print(f'\ndiagonal (mismo residuo): min {min(diagonal):.0f}, '
          f'max {max(diagonal):.0f}, promedio {sum(diagonal) / len(diagonal):.2f}')
    print(f'fuera de la diagonal: min {min(off_diagonal):.0f}, '
          f'max {max(off_diagonal):.0f}, '
          f'promedio {sum(off_diagonal) / len(off_diagonal):.2f}')
    positives = sum(1 for value in off_diagonal if value > 0)
    print(f'sustituciones distintas con puntaje positivo (conservativas): '
          f'{positives} de {len(off_diagonal)}')
