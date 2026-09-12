"""Grafico del dot-plot."""

import matplotlib.pyplot as plt


def plot_dotplot(matrix, seq_rows, seq_cols, title='Dot-plot'):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.imshow(matrix, cmap='Greys', origin='upper', aspect='auto')
    ax.set_title(title)

    if len(seq_cols) <= 60:
        ax.set_xticks(range(len(seq_cols)))
        ax.set_xticklabels(list(seq_cols))
    if len(seq_rows) <= 60:
        ax.set_yticks(range(len(seq_rows)))
        ax.set_yticklabels(list(seq_rows))

    fig.tight_layout()
    plt.show()
