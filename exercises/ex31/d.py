"""3.1d) Programe un algoritmo que recorra un dot-plot (de unos y ceros) y
obtenga un buen (o el mejor) camino posible.

El camino arranca en la esquina arriba-izquierda y va hacia el limite derecho o
inferior, idealmente hasta la esquina abajo-derecha. Pasos posibles: D, V, H.

Consejo del TP: piense/discuta como recorreria el dot-plot y que decision
tomaria a cada paso, usando el principio de dividir y conquistar. La decision
esta explicada en el docstring de paths.py: el mejor camino hasta una casilla
se arma con el mejor camino hasta sus tres vecinas.
"""

from dotplot import build_dotplot, filter_dotplot
from paths import find_path, path_score, path_steps, step_counts
from sequences import load_pair

WINDOW = 10
THRESHOLD = 6
GAP_PENALTY = -1.0


def describe(dotplot, path, label):
    counts = step_counts(path)
    print(f'{label}:')
    print(f'  arranca en {path[0]} y termina en {path[-1]}')
    print(f'  {len(path)} casillas, pasos D={counts["D"]} V={counts["V"]} H={counts["H"]}')
    print(f'  puntaje (matches sobre el camino): {path_score(dotplot, path):.0f}')
    print(f'  primeros pasos: {path_steps(path)[:40]}...')


def run():
    s1, s2 = load_pair()
    dotplot = build_dotplot(s1, s2)

    path = find_path(dotplot, gap_penalty=GAP_PENALTY)
    describe(dotplot, path, 'sobre el dot-plot crudo')

    filtered = filter_dotplot(dotplot, WINDOW, THRESHOLD, diagonal=True)
    filtered_path = find_path(filtered, gap_penalty=GAP_PENALTY)
    describe(filtered, filtered_path,
             f'\nsobre el dot-plot filtrado (w = {WINDOW}, umbral = {THRESHOLD})')

    return path
