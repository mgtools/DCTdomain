### HOMOLOGY BENCHMARK RESULTS

This folder contains the results from using each homology detection tool on the benchmark datasets.

On each dataset, each tool has a csv file with a label for each pair (1 if homolog, 0 if not) and a score. These scores are used to generate the ROC curves seen in the paper. You can replicate these graphs by running AUC.py from this directory.

Tools:

    - DCTdomain/global with FU predicted domains
    - DCTdomain/global with RecCut predicted domains
    - BLAST
    - CSBlast
    - FASTA
    - HHSearch
    - phmmer
    - UBLAST
    - USEARCH

Datasets:

    - pfam_max50
    - pfam_nomax50
    - pfam_localpfam_nomax50
    - supfam_nomax50
    - gene3d_nomax50

Below is a guide on replicating results from https://academic.oup.com/bioinformatics/article/32/17/2636/2450749?, which we replicated in this analysis. All parameters are default except when noted. Documentation for each tool which will be linked under each section. High e-values (10e9) are used to ensure that all hits are returned.

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