"""Defines the Database class, which is used to manage the SQLite database.

__author__ = "Ben Iovino"
__date__ = "3/18/24"
"""

import datetime
import faiss
import sqlite3
import os
import numpy as np
from io import BytesIO
from fingerprint import Fingerprint


class Database:
    """Initializes .db file and manages database operations.

    Attributes:
        path (str): Path to database file.
        conn (sqlite3.Connection): Connection to database.
        cur (sqlite3.Cursor): Cursor for database.
    """


    def __init__(self, dbfile: str, fafile: str = None):
        """Initializes database file and cursor. If a fasta file is given, sequences are read and
        added to the database if they are not already present.
        
        Args:
            dbfile (str): Path to database file.
            fafile (str): Path to fasta file.
        """

        if fafile and (fafile.endswith('.fa') or fafile.endswith('.fasta')):
            print(f'Reading file: {fafile}')
            self.path = dbfile
            seqs = self.read_fasta(fafile)
            self.init_db(seqs)

        # Fasta file not necessary, open database if it exists
        else:
            if not os.path.exists(dbfile):
                raise FileNotFoundError(f'Database file not found: {dbfile}')
            else:
                print(f'Opening database: {dbfile}')
                self.path = os.path.splitext(dbfile)[0]
                self.conn = sqlite3.connect(f'{self.path}.db')
                self.cur = self.conn.cursor()

    
    def close(self):
        """Closes the database connection.
        """

        print(f'Closing database: {self.path}\n')
        self.conn.close()


    def read_fasta(self, fafile: str) -> dict:
        """Returns a dictionary of sequences from a fasta file. They are sorted by lenth to
        optimize fingerprinting.

        Args:
            fafile (str): Path to fasta file.

        Returns:
            dict: Dictionary of sequences.
        """

        seqs = {}
        with open(fafile, 'r', encoding='utf8') as f:
            for line in f:
                if line.startswith('>'):
                    pid = line.strip().split()[0][1:]
                    seqs[pid] = ''
                else:
                    seqs[pid] += line.strip()
        seqs = {k: v for k, v in sorted(seqs.items(), key=lambda item: len(item[1]))}
        
        return seqs
    

    def init_db(self, seqs: dict):
        """Creates a new database file and fills table with sequences.

        Sequences table is primarily used for taking protein sequences to create the fingerprints
        table. Fingerprints are stored in a separate table to allow ease of access when searching
        using FAISS which returns the indices of the fingerprints in a flat numpy array.

        Args:
            seqs (dict): Dictionary of sequences.
        """

        self.path = os.path.splitext(self.path)[0]
        self.conn = sqlite3.connect(f'{self.path}.db')
        self.cur = self.conn.cursor()

        # Create table for protein sequences
        table = """CREATE TABLE IF NOT EXISTS sequences (
                pid text PRIMARY KEY,
                sequence text NOT NULL,
                length integer NOT NULL,
                fpcount integer NOT NULL
                ); """
        self.cur.execute(table)

        # Create table for domain fingerprints
        table = """CREATE TABLE IF NOT EXISTS fingerprints (
                vid integer PRIMARY KEY,
                domain text NOT NULL,
                fingerprint blob NOT NULL,
                pid text NOT NULL,
                FOREIGN KEY(pid) REFERENCES sequences(pid)
                ); """
        self.cur.execute(table)

        # Create table for metadata
        table = """CREATE TABLE IF NOT EXISTS metadata (
                datetime text PRIMARY KEY,
                seq_num integer NOT NULL,
                avg_len real NOT NULL,
                fp_num integer NOT NULL,
                seqs_fp string NOT NULL
                ); """
        self.cur.execute(table)

        # Insert sequences
        insert = """ INSERT OR IGNORE INTO sequences(pid, sequence, length, fpcount)
            VALUES(?, ?, ?, ?) """
        for pid, seq in seqs.items():
            self.cur.execute(insert, (pid, seq, len(seq), 0))
        self.conn.commit()


    def yield_seqs(self, maxlen: int, cpu: int, dim1: int = 3, dim2: int = 80):
        """Yields sequences from the database.

        Args:
            maxlen (int): Maximum length of total sequence to yield
            cpu (int): Number of cpu cores to use (consequently, hard maximum of seqs to yield)
            dim1 (int): First dimension of quantization.
            dim2 (int): Second dimension of quantization.

        Yields:
            list: List of tuples where elements are protein ID and sequence.
        """

        seqs, curr_len, min_size = [], 0, dim1*dim2
        select = """ SELECT pid, sequence, length FROM sequences WHERE fpcount = 0 """
        rows = self.cur.execute(select).fetchall()
        for row in rows:
            pid, seq, length = row
            if (length-2) * dim2 < min_size:  # Short sequences can't be quantized
                continue
            curr_len += length

            # If list is too large (length/number of seqs), yield and reset
            if (len(seqs) > 1 and curr_len > maxlen) or len(seqs) > cpu:
                last_seq = seqs.pop()
                yield seqs
                curr_len = len(last_seq[1])
                seqs = [(last_seq[0], last_seq[1])]
            seqs.append((pid, seq))  # New sequence

        # Last batch in file may be too large
        if len(seqs) > 1 and curr_len > maxlen:
            last_seq = seqs.pop()
            yield seqs

        # Yield unless there are no sequences selected from database
        try:
            yield [(last_seq[0], last_seq[1])]
        except UnboundLocalError:
            if seqs:  # Only one sequence in fasta/database
                yield seqs
            else:
                print('No sequences to fingerprint!\n')

    
    def get_last_vid(self) -> int:
        """Returns the last fingerprint ID in the database.

        Returns:
            int: Last fingerprint ID.
        """

        select = """ SELECT vid FROM fingerprints ORDER BY vid DESC LIMIT 1 """
        try:
            vid = self.cur.execute(select).fetchone()[0] + 1
        except TypeError:
            vid = 1

        return vid


    def add_fprint(self, fp: Fingerprint, lock, counter):
        """Adds domains and fingerprints to database.

        Args:
            fp (Fingerprint): Fingerprint object to add to database.
            lock (multiprocessing.Lock): Lock for multiprocessing.
            counter (multiprocessing.Value): Value object for unique fingerprint ID's (vid)
        """

        # Convert quantizations to bytes for db storage
        quants = np.array([fp.quants[dom] for dom in fp.domains], dtype=np.int8)
        quants_bytes = BytesIO()
        np.save(quants_bytes, quants, allow_pickle=True)

        # Update sequences table with number of fingerprints to keep track of progress
        update = """ UPDATE sequences SET fpcount = ? WHERE pid = ? """
        self.cur.execute(update, (len(fp.domains), fp.pid))

        # Add each domain and it's fingerprint to the fingerprints table
        insert = """ INSERT INTO fingerprints(vid, domain, fingerprint, pid)
            VALUES(?, ?, ?, ?) """
        for dom, quant in zip(fp.domains, quants):
            quants_bytes = BytesIO()
            np.save(quants_bytes, quant, allow_pickle=True)
            with lock:
                counter.value += 1  # ensure unique vid
                self.cur.execute(insert, (counter.value, dom, quants_bytes.getvalue(), fp.pid))
        self.conn.commit()


    def create_index(self):
        """Creates index of fingerprints for fast querying with FAISS.
        """

        # Load fingerprints
        select = """ SELECT fingerprint FROM fingerprints """
        self.cur.execute(select)
        fps = []
        for row in self.cur:
            fprint = np.load(BytesIO(row[0]), allow_pickle=True)
            fps.append(fprint)

        # Create index
        fps = np.array(fps, dtype=np.int8)
        index = faiss.IndexFlatL2(fps.shape[1])
        index.add(fps)
        faiss.write_index(index, f'{self.path}.index')  # db.path is path w/o extension


    def load_fprints(self, pid: str = '') -> list:
        """Returns a list of fingerprints for a given sequence ID.

        Args:
            pid (str): Protein ID of sequence in database for querying.

        Returns:
            list: List of numpy arrays
        """

        select = """ SELECT vid, fingerprint FROM fingerprints WHERE pid = ? """
        self.cur.execute(select, (pid,))

        # Load fingerprints one at a time to save memory
        fprints = []
        for row in self.cur:
            fprint = np.load(BytesIO(row[1]), allow_pickle=True)
            fprints.append((row[0], fprint))
        
        return fprints
    

    def rename_vid(self):
        """Renames fingerprint ID column in database to be sequential. This is necessary for the
        FAISS index to work properly. ID's can be incorrect if multiprocessing is used to add
        fingerprints to the database.
        """

        # Get all fingerprints
        select = """ SELECT vid FROM fingerprints """
        vids = self.cur.execute(select).fetchall()

        # Rename fingerprints
        update = """ UPDATE fingerprints SET vid = ? WHERE vid = ? """
        for i, vid in enumerate(vids):
            self.cur.execute(update, (i+1, vid[0]))
        self.conn.commit()


    def update_metadata(self):
        """Updates metadata in the database.
        """

        print('Updating metadata...')
        select = """ SELECT COUNT(*) FROM sequences """
        num_seqs = self.cur.execute(select).fetchone()[0]
        select = """ SELECT AVG(length) FROM sequences """  
        avg_len = self.cur.execute(select).fetchone()[0]
        select = """ SELECT SUM(fpcount), COUNT(*) FROM sequences WHERE fpcount > 0 """
        nom_dom, dom_seqs = self.cur.execute(select).fetchone()
        if not nom_dom:
            nom_dom = 0

        # Insert new metadata
        insert = """ INSERT INTO metadata(datetime, seq_num, avg_len, fp_num, seqs_fp)
            VALUES(?, ?, ?, ?, ?) """
        date = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        self.cur.execute(insert, (date, num_seqs, avg_len, nom_dom, f'{dom_seqs}/{num_seqs}'))
        self.conn.commit()
        self.db_info()


    def db_info(self):
        """Prints information about the database, updating metadata if necessary.
        """

        select = """ SELECT * FROM metadata ORDER BY datetime DESC LIMIT 1 """
        metadata = self.cur.execute(select).fetchone()
        try:
            print(f'Last Updated: {metadata[0]}')
            print(f'Number of Sequences: {metadata[1]}')
            print(f'Average Sequence Length: {metadata[2]:.2f}')
            print(f'Number of Fingerprints: {metadata[3]} ({metadata[4]} fingerprinted)\n')
        except TypeError:
            self.update_metadata()


    def seq_info(self, seq: str):
        """Prints information about a specific sequence.

        Args:
            seq (str): Protein ID of sequence in database.
        """

        # Get sequence from sequences table
        print(f'Protein ID: {seq}')
        select = """ SELECT sequence FROM sequences WHERE pid = ? """
        try:
            sequence = self.cur.execute(select, (seq,)).fetchone()[0]
        except TypeError:
            print('Sequence not found in database\n')
            return
        
        # Get domains from fingerprints table
        select = """ SELECT domain FROM fingerprints WHERE pid = ? """
        domains = self.cur.execute(select, (seq,)).fetchall()
        
        # Print information
        print(f'Sequence: {sequence}')
        if domains:
            print(f'Domains: {", ".join([dom[0] for dom in domains])}\n')
        else:
            print('No domains in database\n')


    def save_fprints(self, file: str):
        """Saves fingerprints to a npz file formatted as pid, idx, dom, fp.

        Args:
            file (str): Path to save file.
        """

        select = """ SELECT pid, domain, fingerprint FROM fingerprints """
        self.cur.execute(select)
        seqs, idxs, doms, fps  = [], [], [], []
        seq, idx = '', 0

        # Only add sequence ID + starting index once, but add each domain and fingerprint
        for i, row in enumerate(self.cur):
            print(i, row)
            fprint = np.load(BytesIO(row[2]), allow_pickle=True)
            if row[0] != seq:
                seq = row[0]
                seqs.append(seq)
                idxs.append(idx)
            doms.append(row[1])
            fps.append(fprint)
            idx += 1

            if i > 50:
                break

        idxs.append(idx)
        np.savez(file, sid=seqs, idx=idxs, dom=doms, dct=fps)
