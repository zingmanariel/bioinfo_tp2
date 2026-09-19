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

    print(f'path of {len(path)} coordinates')
    print(f'first {COORDS_SHOWN}: {path[:COORDS_SHOWN]}')
    print(f'last {COORDS_SHOWN}:  {path[-COORDS_SHOWN:]}')

    plot_dotplot_with_path(dotplot, path,
                           title=f'Path over the raw dot-plot ({len(path)} steps)')

    filtered = filter_dotplot(dotplot, WINDOW, THRESHOLD, diagonal=True)
    plot_dotplot_with_path(filtered, find_path(filtered),
                           title=f'Path over the filtered dot-plot '
                                 f'(w = {WINDOW}, threshold = {THRESHOLD})')
    return path
