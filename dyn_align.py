from dotplot import build_dotplot
from plots import plot_dotplot


class DotPlot:
    """El dot-plot de dos secuencias, con su matriz y su grafico."""

    def __init__(self, seq1, seq2):
        self.seq1 = seq1
        self.seq2 = seq2
        # seq2 va en las filas y seq1 en las columnas, como asume calc_alignment
        self.matrix = build_dotplot(seq2, seq1)

    def plot(self, title='Dot-plot'):
        plot_dotplot(self.matrix, self.seq2, self.seq1, title=title)


GAP = -1


def get_best_move(dot_plot, score_mat, row, col):
    """(puntaje, padre) del mejor paso que llega a (row, col).

    El puntaje de un paso depende de COMO se llega, asi que elegir el padre y
    calcular el puntaje son la misma decision: gana el que maximiza puntaje del
    padre + puntaje del paso.

      DIAG   alinea seq2[row-1] con seq1[col-1] -> lo que diga el dot-plot
             (1 si es match, 0 si es mismatch)
      VERT   letra de seq2 contra un gap        -> GAP
      HORIZ  letra de seq1 contra un gap        -> GAP

    max() se queda con la primera de las empatadas y la diagonal va primera,
    asi que ante un empate gana el paso diagonal.
    """
    match_score = get_cell_score(dot_plot, row, col)
    options = [
        (score_mat[row - 1][col - 1] + match_score, (row - 1, col - 1)),   # DIAG
        (score_mat[row - 1][col] + GAP, (row - 1, col)),                   # VERT
        (score_mat[row][col - 1] + GAP, (row, col - 1)),                   # HORIZ
    ]
    return max(options, key=lambda option: option[0])


def get_cell_score(dot_plot, row, col):
    """Valor del dot-plot para la celda (row-1, col-1); 0 si cae afuera."""
    if row - 1 < len(dot_plot.matrix) and col - 1 < len(dot_plot.matrix[0]):
        return dot_plot.matrix[row - 1][col - 1]
    return 0


def get_dir(path, i):
    """Direccion del paso que llega a path[i] desde su padre.

    path viene al reves (del final al principio), asi que el padre de path[i]
    es path[i + 1].
    """
    if i + 1 >= len(path):
        return "DIAG"
    row, col = path[i]
    parent_row, parent_col = path[i + 1]
    if parent_row == row - 1 and parent_col == col - 1:
        return "DIAG"
    if parent_row == row - 1:
        return "VERT"
    return "HORIZ"


def calc_alignment(seq1, seq2, dot_plot):

    rows = len(dot_plot.matrix)
    cols = len(dot_plot.matrix[0])

    score_mat = [[0] * (cols + 1) for _ in range(rows + 1)]
    best_parent_mat = [[None] * (cols + 1) for _ in range(rows + 1)]

    for row in range(len(score_mat)):
        for col in range(len(score_mat[0])):
            if col == 0 or row == 0:
                score_mat[row][col] = 0
                if col > 0:
                    best_parent_mat[row][col] = (row, col - 1)
                elif row > 0:
                    best_parent_mat[row][col] = (row - 1, col)
            else:
                # el puntaje del paso depende de la direccion, asi que el padre
                # y el puntaje salen de la misma comparacion
                score_mat[row][col], best_parent_mat[row][col] = get_best_move(dot_plot, score_mat, row, col)  # row-1, col-1, if pos is outside of dotplot score is 0

    #----------
    best_path = []

    row = len(score_mat) - 1
    col = len(score_mat[0]) - 1

    while not(row == 0 and col == 0):
        best_path.append((row, col))
        row, col = best_parent_mat[row][col]

    #-----------
    seq1_aligned = []
    seq2_aligned = []


    for i in range(len(best_path)):
        last_dir = get_dir(best_path, i)
        row, col = best_path[i]

        if last_dir == "DIAG":
            seq1_aligned.insert(0, seq1[col-1])
            seq2_aligned.insert(0, seq2[row-1])

        elif last_dir == "VERT":
            seq1_aligned.insert(0, "-")
            seq2_aligned.insert(0, seq2[row-1])
        else:  # last_dir == "HORIZ"
            seq1_aligned.insert(0, seq1[col-1])
            seq2_aligned.insert(0, "-")


    return ''.join(seq1_aligned), ''.join(seq2_aligned)


def main():
    # seq1 = "ACTGA"
    # seq2 = "ACGTA"

    seq1 = "CAGATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTT"
    seq2 = "CAGATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTT"

    dot_plot = DotPlot(seq1, seq2)
    dot_plot.plot()

    alignment = calc_alignment(seq1, seq2, dot_plot)
    print("Alignment:")
    print(alignment[0])
    print(alignment[1])


if __name__ == "__main__":
    main()
