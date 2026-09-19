"""3.2b) Programe un modulo que convierta la matriz de comparacion en un
dot-plot (de 1s y 0s) tomando como parametros: a) un windowsize y b) un umbral
de puntaje para decir que convierte en 0 y que convierte en 1."""

from dotplot import build_comparison_matrix, comparison_to_dotplot, density
from plots import plot_dotplot_grid
from sequences import load_protein_pair
from substitution import load_matrix

WINDOWS = [3, 7, 11]
THRESHOLDS = [5, 12, 20]


def run():
    s1, s2 = load_protein_pair()
    matrix = load_matrix()
    comparison = build_comparison_matrix(s1, s2, matrix)

    print(f'{"":<10}' + ''.join(f'threshold {t:<5}' for t in THRESHOLDS))
    dotplots = {}
    for window in WINDOWS:
        cells = []
        for threshold in THRESHOLDS:
            dotplot = comparison_to_dotplot(comparison, window, threshold)
            cells.append(f'{density(dotplot):.2%}')
            dotplots[f'w = {window}, threshold = {threshold}'] = dotplot
        print(f'w = {window:<5} ' + ''.join(f'{cell:<15}' for cell in cells))

    print('\nThe window sums the scores along its diagonal, so the threshold '
          'reads as "how much accumulated score does a run of w residues '
          'need". A small window + low threshold lets noise through; a large '
          'window + high threshold keeps only the most conserved stretches '
          'and can cut the diagonal into pieces.')

    shown = {label: dotplots[label] for label in dotplots
             if label.startswith(f'w = {WINDOWS[1]}')}
    plot_dotplot_grid(shown, title=f'Protein dot-plot (w = {WINDOWS[1]})')
    return dotplots
