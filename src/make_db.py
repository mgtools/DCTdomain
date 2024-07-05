"""Makes a database of DCT fingerprints from a fasta file of protein sequences.

__author__ = "Ben Iovino"
__date__ = "12/18/23"
"""

import argparse
import datetime
import logging
import os
import torch
import torch.multiprocessing as mp
from multiprocessing import Pool, Lock, Value
from embedding import Model, Batch
from fingerprint import Fingerprint
from database import Database


def queue_cpu(fp: Fingerprint) -> Fingerprint:
    """Predicts domains and quantizes embeddings on cpu. Returns Fingerprint.

    Args:
        queue (mp.Queue): Queue of embeddings to fingerprint

    Returns:
        Fingerprint: Fingerprint object with quantized domains.
    """

    fp.reccut(2.6)
    fp.quantize([3, 80, 3, 80])
    logging.info(f'{datetime.datetime.now()} Fingerprinted {fp.pid}')

    return fp


def fprint_cpu(batch: list, args: argparse.Namespace, db: Database, lock: mp.Lock, counter: mp.Value):
    """Puts batches of Fingerprints in queue to be quantized on cpu(s). Returns a list of
    Fingerprint objects with quantized domains.

    Args:
        batch (list): List of fingerprints to quantize
        args (argparse.Namespace): Command line arguments
        db (Database): Database object connected to SQLite database
        lock (Lock): Lock object for Value counter
        counter (Value): Value object for counting fingerprints
    """

    with Pool(processes=args.cpu) as pool:
        results = pool.starmap(queue_cpu, [(fp,) for fp in batch])
    for fp in results:
        db.add_fprint(fp, lock, counter)


def queue_gpu(rank: int, queue: mp.Queue, args: argparse.Namespace, db: Database, lock: mp.Lock, counter: mp.Value):
    """Moves through queue of sequences to fingerprint and add to database.

    Args:
        rank (int): GPU to load model on
        queue (mp.Queue): queue of batches to embed and transform
        args (argparse.Namespace): explained in main()
        db (Database): Database object connected to SQLite database
        lock (mp.Lock): Lock object for Value counter
        counter (mp.Value): Value object for counting fingerprints
    """

    # Load tokenizer and encoder
    device = torch.device(f'cuda:{rank}')
    model = Model()  # pLM encoder and tokenizer
    model.to_device(device)

    # Embed batches of sequences
    cpu_queue = []
    while True:
        seqs = queue.get()
        if seqs is None:
            break
        batch = Batch(seqs, model, device)
        batch.embed_batch([15, 21], args.maxlen)

        # Create Fingerprint object and add to queue
        for emb in batch.embeds:
            fp = Fingerprint(pid=emb.pid, seq=emb.seq, embed=emb.embed, contacts=emb.contacts)
            cpu_queue.append(fp)

        # If queue is full, start multiprocess fingerprinting
        if len(cpu_queue) >= args.cpu/args.gpu:
            fprint_cpu(cpu_queue, args, db, lock, counter)
            cpu_queue = []

    # Last batch if queue is not full
    if cpu_queue:
        fprint_cpu(cpu_queue, args, db, lock, counter)


def embed_gpu(args: argparse.Namespace, db: Database, lock: mp.Lock, counter: mp.Value):
    """Puts batches of sequences in queue to be embedded on gpu(s).

    Args:
        args (argparse.Namespace): Command line arguments
        db (Database): Database object connected to SQLite database
        lock (Lock): Lock object for Value counter
        counter (Value): Value object for counting fingerprints
    """

    mp_queue = mp.Queue()
    processes = []
    for rank in range(args.gpu):
        proc = mp.Process(target=queue_gpu, args=(rank, mp_queue, args, db, lock, counter))
        proc.start()
        processes.append(proc)
    for seqs in db.yield_seqs(args.maxlen, args.cpu):
        mp_queue.put(seqs)
    for _ in range(args.gpu):  # send None to each process to signal end of queue
        mp_queue.put(None)
    for proc in processes:
        proc.join()
    db.rename_vid()


def embed_cpu(args: argparse.Namespace, db: Database, lock: mp.Lock, counter: mp.Value):
    """Embeds sequences on cpu.
    
    Args:
        args (argparse.Namespace): Command line arguments
        db (Database): Database object connected to SQLite database
        lock (Lock): Lock object for Value counter
        counter (Value): Value object for counting fingerprints
    """

    model = Model()
    device = torch.device('cpu')
    model.to_device(device)
    
    # Embed one sequence at a time (batching on cpu is very slow)
    cpu_queue = []
    for seqs in db.yield_seqs(1, args.cpu):
        batch = Batch(seqs, model, device)
        batch.embed_batch([15, 21], args.maxlen)

        # Create Fingerprint object and add to queue
        for emb in batch.embeds:
            fp = Fingerprint(pid=emb.pid, seq=emb.seq, embed=emb.embed, contacts=emb.contacts)
            cpu_queue.append(fp)

        # If queue is full, start multiprocess fingerprinting
        if len(cpu_queue) >= args.cpu:
            fprint_cpu(cpu_queue, args, db, lock, counter)
            cpu_queue = []

    # Last batch if queue is not full
    if cpu_queue:
        fprint_cpu(cpu_queue, args, db, lock, counter)
    db.rename_vid()


def main():
    """Sequences from a fasta file of protein sequences go through two processes:

    1. Embedding: sequences are embedded using ESM2, recommended to be done on GPU.
    2. Fingerprinting: domains are cut and quantized, performed on CPU.

    The resulting fingerprints are then written to a database file. Increasing the number of
    both available GPUs and CPUs will speed up the process. Make sure --maxlen is set to a
    value that is appropriate for the available memory.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('--fafile', type=str, required=False, help='fasta file to embed')
    parser.add_argument('--dbfile', type=str, required=True, help='db file to write to')
    parser.add_argument('--maxlen', type=int, default=500, help='max sequence length to embed')
    parser.add_argument('--cpu', type=int, default=1, help='number of cpus to use')
    parser.add_argument('--gpu', type=int, required=False, help='number of gpus to use')
    parser.add_argument('--noindex', action='store_true', help='toggle for not creating index')
    parser.add_argument('--nonpz', action='store_true', help='toggle for not saving fingerprints')
    parser.add_argument('--nodom', action='store_true', help='toggle for not writing domains')
    parser.add_argument('--out', type=str, default='', help='print progress to file')
    args = parser.parse_args()

    # Log results to file, otherwise output is quiet
    if args.out:
        logging.basicConfig(level=logging.INFO, filename=args.out,
                             filemode='w', format='%(message)s')
    
    # Get last vid and create lock and counter
    db = Database(args.dbfile, args.fafile)
    vid = db.get_last_vid()
    lock, counter = Lock(), Value('i', vid)

    # Fingerprint sequences
    print('Fingerprinting sequences...\n')
    if args.gpu:
        embed_gpu(args, db, lock, counter)
    else:
        embed_cpu(args, db, lock, counter)
    db.update_metadata()

    # Create index
    if not args.noindex:
        os.environ['OMP_NUM_THREADS'] = str(args.cpu)
        print('Creating index...')
        db.create_index()
    
    # Save fingerprints to npz, doms to txt
    if not args.nonpz:
        db.save_fprints(f'{args.dbfile}-dct.npz')
    if not args.nodom:
        db.save_doms(f'{args.dbfile}.dom')
    db.close()


if __name__ == '__main__':
    main()
