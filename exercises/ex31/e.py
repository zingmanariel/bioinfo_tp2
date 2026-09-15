"""3.1e) Obtenga el camino como una serie de coordenadas y dibujelo sobre el
dot-plot."""

from dotplot import build_dotplot, filter_dotplot
from paths import find_path
from plots import plot_dotplot_with_path
from sequences import load_pair

WINDOW = 10
THRESHOLD = 6
COORDS_SHOWN = 15


def run():
    s1, s2 = load_pair()
    dotplot = build_dotplot(s1, s2)
    path = find_path(dotplot)

    print(f'camino de {len(path)} coordenadas')
    print(f'primeras {COORDS_SHOWN}: {path[:COORDS_SHOWN]}')
    print(f'ultimas {COORDS_SHOWN}:  {path[-COORDS_SHOWN:]}')

    plot_dotplot_with_path(dotplot, path,
                           title=f'Camino sobre el dot-plot crudo ({len(path)} pasos)')

    filtered = filter_dotplot(dotplot, WINDOW, THRESHOLD, diagonal=True)
    plot_dotplot_with_path(filtered, find_path(filtered),
                           title=f'Camino sobre el dot-plot filtrado '
                                 f'(w = {WINDOW}, umbral = {THRESHOLD})')
    return path
