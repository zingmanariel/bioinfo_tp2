"""Entry point for the TP. Runs one exercise at a time, by name:

    python3 main.py 3.1a
    python3 main.py 3.2c
    python3 main.py          # lists the available exercises

The names match the numbering in the TP PDF, so the command you type, the
section in the report and the file on disk all agree.
"""

import sys

from exercises.ex31 import (a as ex31a, b as ex31b, c as ex31c, d as ex31d,
                            e as ex31e, f as ex31f, g as ex31g, h as ex31h,
                            i as ex31i)
from exercises.ex32 import (a as ex32a, b as ex32b, c as ex32c, d as ex32d,
                            e as ex32e, f as ex32f)

EXERCISES = {
    # 3.1 Dot-plot, filtros, camino y alineamiento (ADN)
    "3.1a": ex31a.run,
    "3.1b": ex31b.run,
    "3.1c": ex31c.run,
    "3.1d": ex31d.run,
    "3.1e": ex31e.run,
    "3.1f": ex31f.run,
    "3.1g": ex31g.run,
    "3.1h": ex31h.run,
    "3.1i": ex31i.run,

    # 3.2 Alineamiento de proteinas (matriz de sustitucion)
    "3.2a": ex32a.run,
    "3.2b": ex32b.run,
    "3.2c": ex32c.run,
    "3.2d": ex32d.run,   # opcional
    "3.2e": ex32e.run,   # opcional
    "3.2f": ex32f.run,   # opcional
}


def usage():
    print("uso: python3 main.py <ejercicio>")
    print("ejercicios disponibles:")
    for name in EXERCISES:
        print(f"  {name}")


def main():
    if len(sys.argv) != 2:
        usage()
        return
    name = sys.argv[1]
    if name not in EXERCISES:
        print(f"'{name}' no es un ejercicio conocido.\n")
        usage()
        return
    EXERCISES[name]()


if __name__ == "__main__":
    main()
