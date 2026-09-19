"""3.1c) Programe un filtro que opere sobre el dot-plot obtenido (el window-size
"w" deberia ser una variable a modificar), y reanalice los resultados.

Alternativas a programar y comparar, segun el grafico del TP:
  - ventanas solapadas vs. adyacentes
  - el dot-plot se reemplaza por el de las ventanas (con reduccion) o las
    ventanas se usan solo como filtro (sin reduccion)
  - umbral alto vs. umbral bajo

Y una cuarta, que no esta en el grafico pero cambia todo: contar los matches de
todo el bloque w x w o solo los de su diagonal (ver la nota de dotplot.py).
"""

from dotplot import build_dotplot, density, filter_dotplot, filter_variants
from plots import plot_dotplot_grid
from sequences import load_pair

WINDOW = 10                  # the TP's "w", the variable to tweak
BLOCK_THRESHOLDS = [20, 25, 30]      # counting the whole block (max = w * w)
DIAGONAL_THRESHOLDS = [4, 6, 8]      # counting only the diagonal (max = w)


def run():
    s1, s2 = load_pair()
    dotplot = build_dotplot(s1, s2)
    print(f'unfiltered dot-plot: {density(dotplot):.2%} of cells on\n')

    print(f'--- high vs. low threshold (w = {WINDOW}, overlapping windows) ---')
    print('counting the WHOLE w x w block:')
    for threshold in BLOCK_THRESHOLDS:
        filtered = filter_dotplot(dotplot, WINDOW, threshold)
        print(f'  threshold {threshold:>2} (out of {WINDOW * WINDOW}): '
              f'{density(filtered):.2%} on')
    print('counting only the window\'s DIAGONAL:')
    for threshold in DIAGONAL_THRESHOLDS:
        filtered = filter_dotplot(dotplot, WINDOW, threshold, diagonal=True)
        print(f'  threshold {threshold:>2} (out of {WINDOW}): '
              f'{density(filtered):.2%} on')

    print('\nThe block filter barely de-clutters anything: a window sitting on '
          'the real diagonal has w matches + noise, and any random window '
          'already has w*w/4 of noise, so both come out similar. Counting the '
          'diagonal instead compares w against w/4, and the background '
          'disappears.')

    print(f'\n--- the four combinations from the TP\'s diagram '
          f'(w = {WINDOW}, threshold {DIAGONAL_THRESHOLDS[1]} on the diagonal) ---')
    variants = filter_variants(dotplot, WINDOW, DIAGONAL_THRESHOLDS[1], diagonal=True)
    for label, filtered in variants.items():
        print(f'  {label:<28} matrix {filtered.shape[0]:>3} x {filtered.shape[1]:<3} '
              f'-> {density(filtered):.2%} on')

    plot_dotplot_grid(variants,
                      title=f'Filters (w = {WINDOW}, threshold = {DIAGONAL_THRESHOLDS[1]} '
                            f'on the diagonal)')
    return variants
