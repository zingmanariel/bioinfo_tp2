"""Sequence handling shared by every exercise: reading FASTA files, generating
sequences with a controlled identity/length (for the analysis of 3.1f) and a
couple of small helpers used when comparing two sequences."""

import os
import random

from Bio import SeqIO
from Bio.Data.IUPACData import protein_letters

# The 20 standard amino acids, in alphabetical order (taken from Biopython so
# the list is not hand-typed).
AA = list(protein_letters)

NUCLEOTIDES = ['A', 'C', 'G', 'T']

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, 'data')

# Los pares con los que trabajan los ejercicios (ver ncbi.RECORDS):
#   ADN       beta-globina humana vs. de raton (mismo gen, especies distintas)
#   Proteina  beta-globina vs. delta-globina humanas (muy parecidas)
#   Lejana    mioglobina humana, para comparar contra una proteina lejana
DNA_PAIR = ('seq1.fasta', 'seq2.fasta')
PROTEIN_PAIR = ('prot1.fasta', 'prot2.fasta')
DISTANT_PROTEIN = 'prot3.fasta'


def data_path(filename):
    return os.path.join(DATA_DIR, filename)


# ---------------------------------------------------------------------------
# Reading sequences from disk
# ---------------------------------------------------------------------------

def load_fasta(path):
    """First record of a FASTA file, as an upper case string.

    The TP works with plain strings ("convertirlas a formato String"), so the
    SeqRecord is unwrapped here and never leaves this module.

    If the file is missing, the sequences are re-downloaded from NCBI (they
    ship with the repo, so this is only a fallback).
    """
    if not os.path.exists(path):
        from ncbi import download_missing   # local import: ncbi imports this module
        download_missing()
    record = next(SeqIO.parse(path, 'fasta'))
    return str(record.seq).upper()


def load_fasta_records(path):
    """Every record of a multi-FASTA, as a list of (id, sequence) tuples."""
    return [(r.id, str(r.seq).upper()) for r in SeqIO.parse(path, 'fasta')]


def load_pair(name1=DNA_PAIR[0], name2=DNA_PAIR[1]):
    """The two DNA sequences to compare, read from data/.

    Defaults to the pair that ships with the repo; pass other file names to use
    your own sequences.
    """
    return load_fasta(data_path(name1)), load_fasta(data_path(name2))


def load_protein_pair(name1=PROTEIN_PAIR[0], name2=PROTEIN_PAIR[1]):
    """The two protein sequences of 3.2, read from data/."""
    return load_fasta(data_path(name1)), load_fasta(data_path(name2))


def write_fasta(path, sequence, description='seq'):
    """Write one sequence as FASTA (used to hand sequences over to BLAST)."""
    with open(path, 'w') as f:
        f.write(f'>{description}\n')
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i + 60] + '\n')
    return path


# ---------------------------------------------------------------------------
# Synthetic sequences, to control identity and length (3.1f)
# ---------------------------------------------------------------------------

def generate_random_nt_sequence(length, gc=0.5):
    """Random DNA sequence with the requested GC content (gc=0.5 => uniform)."""
    weights = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]   # A, C, G, T
    return ''.join(random.choices(NUCLEOTIDES, weights=weights, k=length))


def generate_random_aa_sequence(length):
    """Random protein sequence, every residue equally likely."""
    return ''.join(random.choices(AA, k=length))


def mutate(sequence, identity, alphabet=None):
    """Copy of `sequence` where a fraction `identity` of the positions is kept.

    Substitutions only (no indels), so both sequences keep the same length.
    Handy for 3.1f: build pairs at 95%, 70%, 40% identity and look at how the
    dot-plot and the path degrade.
    """
    alphabet = alphabet or NUCLEOTIDES
    mutated = []
    for symbol in sequence:
        if random.random() < identity:
            mutated.append(symbol)
        else:
            others = [s for s in alphabet if s != symbol]
            mutated.append(random.choice(others))
    return ''.join(mutated)


def insert_indels(sequence, n, max_size=5, alphabet=None):
    """Apply `n` random insertions/deletions, so the pair also differs in
    length (the case that forces the path off the main diagonal)."""
    alphabet = alphabet or NUCLEOTIDES
    result = sequence
    for _ in range(n):
        size = random.randint(1, max_size)
        position = random.randint(0, max(0, len(result) - size))
        if random.random() < 0.5:
            result = result[:position] + result[position + size:]
        else:
            insertion = ''.join(random.choices(alphabet, k=size))
            result = result[:position] + insertion + result[position:]
    return result


# ---------------------------------------------------------------------------
# Comparing two sequences
# ---------------------------------------------------------------------------

def identity(row1, row2):
    """Fraction of identical positions between two rows of equal length.

    Meant for aligned rows (gaps included); a gap never counts as a match.
    """
    if len(row1) != len(row2):
        raise ValueError('Both rows must have the same length.')
    matches = sum(1 for a, b in zip(row1, row2) if a == b and a != '-')
    return matches / len(row1) if row1 else 0.0
