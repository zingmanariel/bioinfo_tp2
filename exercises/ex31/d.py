"""3.1d) Programe un algoritmo que recorra un dot-plot (de unos y ceros) y
obtenga un buen (o el mejor) camino posible.

El camino arranca en la esquina arriba-izquierda y va hacia el limite derecho o
inferior, idealmente hasta la esquina abajo-derecha. Pasos posibles: D, V, H.

Consejo del TP: piense/discuta como recorreria el dot-plot y que decision
tomaria a cada paso, usando el principio de dividir y conquistar.
"""

from dotplot import build_dotplot, filter_dotplot
from paths import find_path, path_score, path_steps, step_counts
from sequences import load_pair

WINDOW = 10
THRESHOLD = 6


def describe(dotplot, path, label):
    counts = step_counts(path)
    print(f'{label}:')
    print(f'  starts at {path[0]} and ends at {path[-1]}')
    print(f'  {len(path)} cells, steps D={counts["D"]} V={counts["V"]} H={counts["H"]}')
    print(f'  score (matches along the path): {path_score(dotplot, path):.0f}')
    print(f'  first steps: {path_steps(path)[:40]}...')


def run():
    s1, s2 = load_pair()
    dotplot = build_dotplot(s1, s2)

    path = find_path(dotplot)
    describe(dotplot, path, 'over the raw dot-plot')

    filtered = filter_dotplot(dotplot, WINDOW, THRESHOLD, diagonal=True)
    filtered_path = find_path(filtered)
    describe(filtered, filtered_path,
             f'\nover the filtered dot-plot (w = {WINDOW}, threshold = {THRESHOLD})')

    return path
