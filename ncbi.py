"""Everything that talks to NCBI: downloads the sequences the exercises work
with and caches them in data/, so they don't need to be requested again.

The sequences already ship downloaded in the repo. If you delete data/, the
first run downloads them again.
"""

import os

from Bio import Entrez, SeqIO

from sequences import DATA_DIR, data_path

# NCBI asks for a contact address on every request.
Entrez.email = 'arielz@gmail.com'

# Which sequence goes in each file under data/.
#
# DNA: the same gene (beta-globin) in human and in mouse. They're homologs,
# so the dot-plot shows a clear diagonal with divergent regions, which is
# exactly what 3.1 needs.
#
# Protein: human beta-globin against human delta-globin (very similar) and
# against myoglobin (distantly similar), for the analysis in 3.2c.
RECORDS = {
    'seq1.fasta': ('NM_000518.5', 'nucleotide', 'human HBB (mRNA)'),
    'seq2.fasta': ('NM_008220.3', 'nucleotide', 'mouse Hbb-bs (mRNA)'),
    'prot1.fasta': ('NP_000509.1', 'protein', 'human beta-globin'),
    'prot2.fasta': ('NP_000510.1', 'protein', 'human delta-globin'),
    'prot3.fasta': ('NP_005359.1', 'protein', 'human myoglobin'),
}


def fetch(accession, db, filename):
    """Downloads one record from NCBI and saves it as FASTA in data/filename."""
    print(f'Downloading {accession} from NCBI...')
    handle = Entrez.efetch(db=db, id=accession, rettype='fasta', retmode='text')
    record = next(SeqIO.parse(handle, 'fasta'))
    handle.close()

    os.makedirs(DATA_DIR, exist_ok=True)
    path = data_path(filename)
    SeqIO.write(record, path, 'fasta')
    return path


def download_missing():
    """Downloads the files under data/ that are missing, and returns which ones."""
    downloaded = []
    for filename, (accession, db, description) in RECORDS.items():
        if os.path.exists(data_path(filename)):
            continue
        fetch(accession, db, filename)
        downloaded.append(f'{filename} ({description})')
    return downloaded


if __name__ == '__main__':
    missing = download_missing()
    print(f'Downloaded: {missing}' if missing else 'data/ was already complete.')
