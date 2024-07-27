"""Queries CATH20 database with sequences from cath20_queries.txt using DCTdomain.

__author__ = "Ben Iovino"
__date__ = "4/16/24"
"""

import argparse
import faiss
import logging
import os
import sys
sys.path.append(os.getcwd()+'/src')  # Add src to path
import numpy as np
from src.database import Database
from src.query_db import get_top_hits


def get_queries(file: str) -> dict[str, str]:
    """Returns a dictionary of query PID's and their classification.

    Args:
        file (str): Path to queries file.
    
    Returns:
        dict[str, str]: key: PID, value: classification
    """

    queries = {}
    with open(file, 'r', encoding='utf8') as file:
        for line in file:
            if line.startswith('>'):
                line = line.split('|')
                queries[line[0][1:]] = line[1].strip()  # key: PID, value: classification

    return queries


def search_db(path: str, queries: dict[str, str], khits: int, metric):
    """Similar to search_db in query_db.py, but queries benchmark db with sequences from
    queries text file (which already exist in the db).

    Args:
        path (str): Path to database and index files
        queries (dict[str, str]): Dictionary of query PID's and their homologous superfamily
        classification
        khits (int): Number of hits to return
        metric (str): Distance metric used in faiss index
    """

    db = Database(path)
    index = faiss.read_index(path.replace('.db', '.index'))
    if metric == 'l1':
        index.metric_type = faiss.METRIC_L1
    if metric == 'ip':
        index.metric_type = faiss.METRIC_INNER_PRODUCT
    for pid, fam in queries.items():
        qfps = db.load_fprints(pid=f'{pid}|{fam}')
        que_arrs = np.array([fp[1] for fp in qfps])
        que_ind = np.array([fp[0] for fp in qfps])
        dm, im = index.search(que_arrs, khits)
        get_top_hits(dm, im, khits, db, db, que_ind, metric)

    db.close()


def main():
    """
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('--cpu', type=int, default=1, help='number of cpus to use for search')
    parser.add_argument('--khits', type=int, default=300, help='number of hits to return')
    args = parser.parse_args()

    path = 'bench/cathdb/data'
    query = 'cath20_queries.fa'
    db = 'cath20.db'

    # Read queries sequences
    queries = get_queries(f'{path}/{query}')
    logging.basicConfig(level=logging.INFO, filename=f'{path}/results_dct.txt',
                         filemode='w', format='%(message)s')
    
    # Get queries from CATH20 database and search against database
    os.environ['OMP_NUM_THREADS'] = str(args.cpu)
    search_db(f'{path}/{db}', queries, args.khits, 'l1')
        

if __name__ == '__main__':
    main()
