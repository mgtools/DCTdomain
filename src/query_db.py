"""Queries a database of DCT fingerprints for most similar fingerprints to each query sequence.

__author__ = "Ben Iovino"
__date__ = "3/19/23"
"""

import argparse
import faiss
import logging
import os
import multiprocessing as mp
import numpy as np
from database import Database
from make_db import embed_cpu, embed_gpu


def get_top_hits(dm: np.ndarray, im: np.ndarray, top: int, fp_db: Database, query_db: Database, que_ind: list, metric):
    """Logs the top hits for each query sequence. Distance and index matrices are of dimensions
    (number of query fingerprints x khits). Indices correspond to vector ID's (vid) in the
    fingerprint databases.

    Args:
        dm (np.ndarray): Distance matrix
        im (np.ndarray): Index matrix
        top (int): Number of top hits to return
        fp_db (Database): Database object for fingerprint database
        query_db (Database): Database object for query database
        que_ind (list): List of query vector ID's, allows for querying db against itself
        metric (str): Distance metric used in faiss index
    """

    # Store all hits in a dict and sort by distance
    top_hits = {}
    for i, khits in enumerate(dm):
        for j, dist in enumerate(khits):
            top_hits[i, j] = dist

    # Sort by distance or similarity
    if metric == ('l1' or 'l2'):
        top_hits = dict(sorted(top_hits.items(), key=lambda x: x[1]))
    if metric == 'ip':  # this function used in cathdb benchmark, sort with sim values
        top_hits = dict(sorted(top_hits.items(), key=lambda x: -1*x[1]))

    # Get protein ID and domain for each hit (query and db)
    for i, index in enumerate(list(top_hits.keys())[:top]):

        # Database and query vector ID's (0-indexed in faiss, 1-indexed in SQLite db)
        db_vid = int(im[index[0], index[1]])
        if db_vid == -1:  # No more hits
            break
        q_vid = int(que_ind[index[0]])

        # Log results
        select = """ SELECT pid, domain FROM fingerprints WHERE vid = ? """
        db_pid, db_domain = fp_db.cur.execute(select, (db_vid+1,)).fetchone()
        q_pid, q_domain = query_db.cur.execute(select, (q_vid,)).fetchone()  # not taken from faiss
        logging.info('Query: %s %s, Result %s: %s %s, Similarity: %s',
                      q_pid, q_domain, i+1, db_pid, db_domain, 1-(top_hits[index]/17000))


def search_db(args: argparse.Namespace, query_db: str, fp_db: str, metric: str = 'l1'):
    """Searches a database of DCT fingerprints for the most similar protein to each query sequence.

    Args:
        args (argparse.Namespace): Command line arguments
        query_db (str): Name of query database
        fp_db (str): Name of fingerprint database
        metric (str): Distance metric used in faiss index (default 'l1')
    """

    # Connect to databases
    query_db = Database(query_db)
    print('Loading index...\n')
    index = faiss.read_index(fp_db.replace('.db', '.index'))
    index.metric_type = faiss.METRIC_L1
    fp_db = Database(fp_db)

    # Get each sequence from query db and compare to db
    select = """ SELECT pid FROM sequences """
    query_fps = query_db.cur.execute(select).fetchall()
    print('Querying database...\n')
    for query in query_fps:
        qfps = query_db.load_fprints(pid=query[0])
        que_arr = np.array([fp[1] for fp in qfps])
        que_ind = np.array([fp[0] for fp in qfps])
        dm, im = index.search(que_arr, args.khits)  # distance, index matrices
        get_top_hits(dm, im, args.khits, fp_db, query_db, que_ind, metric)

    query_db.close()
    fp_db.close()


def main():
    """Processes sequences same as make_db.py and queries --db for top --khits fingerprints for
    each sequence in the query database.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('--query', type=str, required=True, help='can be .fa or .db file')
    parser.add_argument('--db', type=str, required=True, help='fingerprint database (.db)')
    parser.add_argument('--out', type=str, default=False, help='output file')
    parser.add_argument('--maxlen', type=int, default=500, help='max sequence length to embed')
    parser.add_argument('--khits', type=int, default=100, help='number of hits to return')
    parser.add_argument('--cpu', type=int, default=1, help='number of cpus to use')
    parser.add_argument('--gpu', type=int, default=False, help='number of gpus to use')
    args = parser.parse_args()

    # Logging for either stdout or file
    if args.out:
        logging.basicConfig(level=logging.INFO, filename=args.out,
                             filemode='w', format='%(message)s')
    else:
        logging.basicConfig(level=logging.INFO, format='%(message)s')

    # Embed query sequences
    query_db = os.path.splitext(args.query)[0]+'.db'
    db = Database(query_db, args.query)
    vid = db.get_last_vid()
    lock, counter = mp.Lock(), mp.Value('i', vid)
    if args.gpu:
        embed_gpu(args, db, lock, counter)
    else:
        embed_cpu(args, db, lock, counter)
    db.update_metadata()

    # Query database for most similar sequence
    os.environ['OMP_NUM_THREADS'] = str(args.cpu)
    search_db(args, query_db, args.db)
   

if __name__ == '__main__':
    main()
