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

    print(f'{"":<10}' + ''.join(f'umbral {t:<8}' for t in THRESHOLDS))
    dotplots = {}
    for window in WINDOWS:
        cells = []
        for threshold in THRESHOLDS:
            dotplot = comparison_to_dotplot(comparison, window, threshold)
            cells.append(f'{density(dotplot):.2%}')
            dotplots[f'w = {window}, umbral = {threshold}'] = dotplot
        print(f'w = {window:<5} ' + ''.join(f'{cell:<15}' for cell in cells))

    print('\nLa ventana suma los puntajes a lo largo de su diagonal, asi que el '
          'umbral se lee como "cuanto puntaje acumulado exige una corrida de w '
          'residuos". Ventana chica + umbral bajo deja pasar el ruido; ventana '
          'grande + umbral alto se queda solo con los tramos mas conservados y '
          'puede cortar la diagonal en pedazos.')

    shown = {label: dotplots[label] for label in dotplots
             if label.startswith(f'w = {WINDOWS[1]}')}
    plot_dotplot_grid(shown, title=f'Dot-plot de proteinas (w = {WINDOWS[1]})')
    return dotplots
