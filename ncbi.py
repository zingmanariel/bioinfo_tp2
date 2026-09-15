"""Todo lo que habla con NCBI: baja las secuencias con las que trabajan los
ejercicios y las cachea en data/, para no volver a pedirlas.

Las secuencias ya vienen descargadas en el repo. Si borra data/, la primera
corrida las vuelve a bajar.
"""

import os

from Bio import Entrez, SeqIO

from sequences import DATA_DIR, data_path

# NCBI pide una direccion de contacto en cada request.
Entrez.email = 'arielz@gmail.com'

# Que secuencia va en cada archivo de data/.
#
# ADN: el mismo gen (beta-globina) en humano y en raton. Son homologos, asi que
# el dot-plot muestra una diagonal clara con regiones divergentes, que es
# justo lo que hace falta para 3.1.
#
# Proteinas: beta-globina humana contra delta-globina humana (muy parecidas) y
# contra mioglobina (parecidas de lejos), para el analisis de 3.2c.
RECORDS = {
    'seq1.fasta': ('NM_000518.5', 'nucleotide', 'HBB humano (mRNA)'),
    'seq2.fasta': ('NM_008220.3', 'nucleotide', 'Hbb-bs raton (mRNA)'),
    'prot1.fasta': ('NP_000509.1', 'protein', 'beta-globina humana'),
    'prot2.fasta': ('NP_000510.1', 'protein', 'delta-globina humana'),
    'prot3.fasta': ('NP_005359.1', 'protein', 'mioglobina humana'),
}


def fetch(accession, db, filename):
    """Baja un record de NCBI y lo guarda como FASTA en data/filename."""
    print(f'Descargando {accession} de NCBI...')
    handle = Entrez.efetch(db=db, id=accession, rettype='fasta', retmode='text')
    record = next(SeqIO.parse(handle, 'fasta'))
    handle.close()

    os.makedirs(DATA_DIR, exist_ok=True)
    path = data_path(filename)
    SeqIO.write(record, path, 'fasta')
    return path


def download_missing():
    """Baja los archivos de data/ que no esten, y devuelve cuales bajo."""
    downloaded = []
    for filename, (accession, db, description) in RECORDS.items():
        if os.path.exists(data_path(filename)):
            continue
        fetch(accession, db, filename)
        downloaded.append(f'{filename} ({description})')
    return downloaded


if __name__ == '__main__':
    missing = download_missing()
    print(f'Descargados: {missing}' if missing else 'data/ ya estaba completo.')
