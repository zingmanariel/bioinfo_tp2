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

WINDOW = 10                  # el "w" del TP, la variable a modificar
BLOCK_THRESHOLDS = [20, 25, 30]      # contando todo el bloque (max = w * w)
DIAGONAL_THRESHOLDS = [4, 6, 8]      # contando solo la diagonal (max = w)


def run():
    s1, s2 = load_pair()
    dotplot = build_dotplot(s1, s2)
    print(f'dot-plot sin filtrar: {density(dotplot):.2%} de casillas prendidas\n')

    print(f'--- umbral alto vs. bajo (w = {WINDOW}, ventanas solapadas) ---')
    print('contando TODO el bloque w x w:')
    for threshold in BLOCK_THRESHOLDS:
        filtered = filter_dotplot(dotplot, WINDOW, threshold)
        print(f'  umbral {threshold:>2} (de {WINDOW * WINDOW}): '
              f'{density(filtered):.2%} prendidas')
    print('contando solo la DIAGONAL de la ventana:')
    for threshold in DIAGONAL_THRESHOLDS:
        filtered = filter_dotplot(dotplot, WINDOW, threshold, diagonal=True)
        print(f'  umbral {threshold:>2} (de {WINDOW}): '
              f'{density(filtered):.2%} prendidas')

    print('\nEl filtro de bloque casi no descongestiona: una ventana sobre la '
          'diagonal real tiene w matches + ruido, y una ventana cualquiera ya '
          'tiene w*w/4 de ruido, asi que ambas dan parecido. Contando la '
          'diagonal, en cambio, se comparan w contra w/4 y el fondo desaparece.')

    print(f'\n--- las cuatro combinaciones del grafico del TP '
          f'(w = {WINDOW}, umbral {DIAGONAL_THRESHOLDS[1]} en la diagonal) ---')
    variants = filter_variants(dotplot, WINDOW, DIAGONAL_THRESHOLDS[1], diagonal=True)
    for label, filtered in variants.items():
        print(f'  {label:<28} matriz {filtered.shape[0]:>3} x {filtered.shape[1]:<3} '
              f'-> {density(filtered):.2%} prendidas')

    plot_dotplot_grid(variants,
                      title=f'Filtros (w = {WINDOW}, umbral = {DIAGONAL_THRESHOLDS[1]} '
                            f'sobre la diagonal)')
    return variants
