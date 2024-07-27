The bench/ directory contains code to reproduce the results and figures reported in the paper. Each subdirectory contains a readme for how to run the benchmark.

To run these programs, install and activate the necessary conda environment as such:

```
cd bench/
conda env create --file env_bench.yml
conda activate dctdomain_bench
```

To use USEARCH, you can download the binary from the website (https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz) and add it to your PATH.

## HOMOLOGY DETECTION BENCHMARKS
Below is a guide on how to use each tool we benchmarked for Figures 4 and 5. All parameters are default except when noted. Documentation for each tool which will be linked under each section. High e-values (10e9) are used to ensure that all hits are returned. For each tool, we performed a pairwise comparison between each sequence in each dataset. The results are stored in the homo/results/ directory.

### BLAST
https://www.ncbi.nlm.nih.gov/books/NBK52636/

To make and search a database:
```
makeblastdb -in <db.fasta> -dbtype <prot/nucl> -parse_seqids -out <db_name>

blastp -query <query.fasta> -db <db_name> -evalue 1000000000
```

### CSBLAST (REQUIRES LEGACY BLAST)
https://github.com/cangermueller/csblast

K4000.lib and K4000.crf files come with csblast installation. The K4000.lib
file is used in this project. --blast-path should be set to wherever
blast-legacy is installed. 

If you are installing csblast and legacy-blast with conda, the K4000 files will be in the ~./anaconda3/data directory and blast-legacy executables will be in the ~./anaconda3/envs/\<env>\/bin directory.

```
formatdb -t <db> -i <db.fasta> -p T -l formatdb.log

csblast -i <query.fasta> -d <db.fasta> -D <K4000.lib/crf> --blast-path <blastpgp dir> -e 1000000000
```

### FASTA
https://fasta.bioch.virginia.edu/wrpearson/fasta/fasta_guide.pdf

Search a query against another sequence:
```
fasta36 <query.fasta> <db.fasta> -b 1000000000'
```

### HHSEARCH
https://github.com/soedinglab/hh-suite

If you are going to utilize profiles, first make sure you have a database downloaded from one of these sources:
1) https://uniclust.mmseqs.com/
2) https://bfd.mmseqs.com/
3) https://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/

This project uses the scop40 database from the third link because it is small. The original paper performing this analysis used 'hhmake' instead, but this decreases the sensitivity of the search.

To make a database, run the following commands:

```
ffindex_from_fasta -s <db>_fas.ffdata <db>_fas.ffindex <db>.fas

hhblits_omp -i <db>_fas -d <scop40/bfd/uniclust30> -oa3m <db>_a3m -n 2 -cpu 1 -v 0

ffindex_apply <db>_a3m.ffdata <db>_a3m.ffindex -i <db>_hmm.ffindex -d <db>_hmm.ffdata -- hhmake -i stdin -o stdout -v 0

cstranslate -f -x 0.3 -c 4 -I a3m -i <db>_a3m -o <db>_cs219
```

To search against this database:

```
hhsearch -i <query.fasta> -d <db> -E 1000000000'
```

### PHMMER
http://eddylab.org/software/hmmer/Userguide.pdf (phmmer instructions on p 41/227)

phmmer reports negative bit scores, so we set the bit score of every search that reported in no hits to be -200.

To search a sequence against another:

```
phmmer --max -E 1000000000 <query.fasta> <db.fasta>
```

### UBLAST
https://www.drive5.com/usearch/manual/cmd_ublast.html
https://drive5.com/usearch/manual/output_files.html

Both ublast and usearch can have different output formats, but we only considered the bit score in this analysis. More options can be found in the second link above.

To make and search a database:
```
usearch -makeudb_ublast <db.fasta> -output <db.udb>

usearch -ublast <query.fasta> -db <db.fasta> -evalue 1e-9 -userout hits.txt -userfields bits
```

### USEARCH
https://drive5.com/usearch/manual/cmd_usearch_local.html

To search a query against another:
```
usearch -search_local <query.fasta> -db <db.fasta> -evalue 1e-9 -userout hits.txt -userfields bits
```

### PROST
https://github.com/MesihK/prost/

prostDistance outputs a positive distance value that is smaller the more similar two sequence are. This is the opposite of bit scores, so we multiply the distance by -1 to make it comparable to the scoring systems of the other tools for AUC calculations.


To search a query against another:
```
db_seq = pyprost.quantSeq(<db.fasta>)

query_seq = pyprost.quantSeq(<query.fasta>)

dist = -1 * pyprost.prostDistance(query_seq, db_seq)
```

## HOMOLOGY SEARCH BENCHMARKS
Below is a guide on how to use each tool we benchmarked for Figure X. All parameters are default except when noted. Documentation for each tool which will be linked under each section. For each tool, we had a set of queries that we searched against a target database.

### MMseqs2
https://github.com/soedinglab/mmseqs2/wiki

To make a database (for both query and target):
```
mmseqs createdb <db>.fa <db>DB
```

And an index for the target database:
```
mmseqs createindex <db>DB -k 7 --split 1
```

To search a query against a target:
```
mmseqs search <query>DB <db>DB <results> tmp --min-ungapped-score 0 -e 10000 -k 7 --k-score 80 --max-seqs 4000 -s 7.5
mmseqs convertalis <query>DB <db>DB <results> <results>.txt
```

### KnnProtT5
We did not use their code, but we generated mean embeddings with ProtT5 and ESM-2 (separately) and used FAISS's flat index and K-NN search algorithm with cosine similarity to find the nearest neighbors for each query sequence. The code for this can be found in cathdb/run_mean.py.
