"""Construccion de la matriz de dot-plot.

Convencion: build_dotplot(seq_rows, seq_cols) arma una matriz de
len(seq_rows) x len(seq_cols), con 1 en cada match y 0 en cada mismatch.
"""


def build_dotplot(seq_rows, seq_cols):
    return [
        [1 if row_char == col_char else 0 for col_char in seq_cols]
        for row_char in seq_rows
    ]
