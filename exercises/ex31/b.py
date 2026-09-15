"""3.1b) Arme un modulo para graficar el dot-plot.

Analice los resultados: para secuencias largas de ADN el dot-plot se ve muy
"congestionado" (ese es el problema que ataca el filtro de 3.1c).
"""

from dotplot import build_dotplot, density
from plots import plot_dotplot
from sequences import load_pair

SMALL = 40


def run():
    s1, s2 = load_pair()

    small = build_dotplot(s1[:SMALL], s2[:SMALL])
    plot_dotplot(small, s1[:SMALL], s2[:SMALL],
                 title=f'Dot-plot de las primeras {SMALL} bases')

    dotplot = build_dotplot(s1, s2)
    plot_dotplot(dotplot, title=f'Dot-plot completo ({len(s1)} x {len(s2)} bases)')

    print(f'{density(dotplot):.1%} de las casillas prendidas: la diagonal '
          f'existe pero queda tapada por el ruido de fondo.')
    return dotplot
