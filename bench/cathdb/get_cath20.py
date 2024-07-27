"""Downloads CATH20 dataset and prepares it for benchmarking.

__author__ = "Ben Iovino"
__date__ = "4/13/24"
"""

import argparse
import os
import sys
sys.path.append(os.getcwd()+'/src')  # Add src to path
import subprocess as sp
from urllib.request import urlretrieve
from zipfile import ZipFile


def download_file(url: str, filename: str, path: str):
    """Downloads filename (any file) from url to path.

    Args:
        url (str): URL to download.
        filename (str): Name of file to save to.
        path (str): Path to save file to.
    """

    if not os.path.exists(path):
        os.makedirs(path)
    print(f'Downloading {url}/{filename} to {path}...')
    urlretrieve(f'{url}/{filename}', f'{path}/{filename}')

    # Unzip if necessary
    if filename.endswith('.gz'):
        with ZipFile(f'{path}/{filename}', 'r') as zip_ref:
            zip_ref.extractall(path)
        os.remove(f'{path}/{filename}')


def read_classes(path: str) -> dict[str, str]:
    """Returns dictionary containing classification for every protein in CATH20.

    Args:
        path (str): Path that contains cath20.fa and cath20_class.txt.

    Returns:
        dict[str, str]: key: CATH ID, value: classification
    """

    classes: dict[str, str] = {}  # key: CATH ID, value: classification
    with open(f'{path}/cath-domain-list-v4_2_0.txt', 'r', encoding='utf8') as file:
        for line in file:
            if line.startswith('#'):
                continue
            line = line.split()
            classes[line[0]] = '.'.join(line[1:5])

    return classes


def modify_fasta(path: str, classes: dict[str, str]):
    """Modifies headers in cath20.fa to include homologous superfamily classification.

    Args:
        path (str): Path that contains cath20.fa.
        classes (dict[str, str]): key: CATH ID, value: classification
    """

    # Read fasta file for CATH IDs and sequences
    seqs: dict[str, str] = {}  # key: "CATH ID|CATH CLASS", value: sequence
    queries: dict[str, list: str] = {}  # key: classification, value: list of pid's in same class
    with open(f'{path}/cath-dataset-nonredundant-S20-v4_2_0.fa', 'r', encoding='utf8') as file:
        for line in file:
            if line.startswith('>'):
                cath_id = line.split('|')[2].split('/')[0]  # i.e. 16vpA00
                class_id = classes[cath_id]  # i.e. 3.30.930.10
                queries[class_id] = queries.get(class_id, []) + [cath_id]
            else:
                seqs[f'{cath_id}|{class_id}'] = line.strip()  # i.e. 16vpA00|3.30.930.10

    # Rewrite fasta file with modified headers
    with open(f'{path}/cath20.fa', 'w', encoding='utf8') as file:
        for pid, seq in seqs.items():
            file.write(f'>{pid}\n{seq}\n')

    # Write queries to file
    with open(f'{path}/cath20_queries.fa', 'w', encoding='utf8') as file:
        for class_id, cath_ids in queries.items():
            if len(cath_ids) > 1:  # Ignore families with only one member
                for cath_id in cath_ids:  # Write each protein in family
                    pid = f'{cath_id}|{class_id}'
                    file.write(f'>{cath_id}|{class_id}\n{seqs[pid]}\n')


def main():
    """Downloads CATH20 v4.2.0 and fingerprints sequences.

    CATH List File (CLF) Format 2.0
    -------------------------------
    Column 1:  CATH domain name (seven characters)
    Column 2:  Class number
    Column 3:  Architecture number
    Column 4:  Topology number
    Column 5:  Homologous superfamily number
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('--maxlen', type=int, default=500, help='max sequence length to embed')
    parser.add_argument('--cpu', type=int, default=1, help='number of cpus to use')
    parser.add_argument('--gpu', type=int, required=False, help='number of gpus to use')
    args = parser.parse_args()

    # Download files and prepare sequences
    path = 'bench/cathdb/data'
    if not os.path.exists(path):
        url = 'ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/all-releases/v4_2_0'
        download_file(f'{url}/non-redundant-data-sets', 'cath-dataset-nonredundant-S20-v4_2_0.fa', path)
        download_file(f'{url}/cath-classification-data', 'cath-domain-list-v4_2_0.txt', path)
        classes = read_classes(path)
        modify_fasta(path, classes)

    # Fingerprint fasta file
    if args.gpu:
        sp.run(['python', 'src/make_db.py', f'--fafile={path}/cath20.fa', f'--dbfile={path}/cath20',
                f'--maxlen={args.maxlen}', f'--cpu={args.cpu}', f'--gpu={args.gpu}'])
    else:
        sp.run(['python', 'src/make_db.py', f'--fafile={path}/cath20.fa', f'--dbfile={path}/cath20',
            f'--maxlen={args.maxlen}', f'--cpu={args.cpu}'])


if __name__ == "__main__":
    main()
