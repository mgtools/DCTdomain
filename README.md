# DCTdomain
Protein domain-level embedding via ESM-2 and DCT for homology detection and similarity search.

Citation: Benjamin G. Iovino, Haixu Tang and Yuzhen Ye. 2024. Protein domain embeddings for fast and accurate similarity search. Genome Research (doi: 10.1101/gr.279127.124).

## Installation
With git and conda, you can clone this repository and install the required dependencies with the following commands:

```
git clone https://github.com/mgtools/DCTdomain
cd DCTdomain
conda env create --file env.yml
conda activate dctdomain
g++ -o src/RecCut src/RecCut.cpp
```
**Make sure that you compile RecCut as above before you start to use DCTdomain.**

## Approach
An overview of the pipeline is visualized in the figure below:

![This is an image](https://github.com/mgtools/DCTdomain/blob/main/misc/DCTdomain-diag.png)

## Generating DCT fingerprints
To create a database of DCT fingerprints for which you can query against, run the following command:

```
python src/make_db.py --fafile <.fa file> --dbfile <output>
```

This will generate a SQLite database file with two tables containing, one for protein sequences and one for DCT fingerprints. If the process is interrupted for any reason, you can restart it and every sequence that has already had fingerprints generated will be skipped.

You can interact with the database like any other SQLite database, but database.py provides a simple interface to interact with the database. For example, the following command will print basic information about the database, such the number of sequences, average length of sequences, and the number of fingerprints in the database:

```
from database import Database
db = Database('<.db file>')
db.db_info()
db.close()
```

You can also search the database for a specific sequence with the following command:

```
db.seq_info('<fasta id>')
```

Alongside the database file, a .npz file will be generated as well for easy access to the DCT fingerprints. The .npz file has four arrays: 'sids' containing the sequence ID's, 'dom' containing the domain predictions for each fingerprint, 'dct' containing the DCT fingerprints corresponding to each domain, and 'idx' containing the index of the first domain/fingerprint in the 'dom' and 'dct' arrays corresponding to each sid. If you do not need the .npz file, you can use the `--nonpz` option to skip generating it.

An .index file will also be generated to allow for faster searching of the database using FAISS. This index is a flat index and allows for KNN search with L1 distance. If you do not need the .index file, you can use the `--noindex` option to skip generating it.

## Multiprocessing and GPU support
Depending on the length of the protein sequences in your fasta file, you may need to adjust the `--maxlen` parameter when calling make_db.py. The default value is 500, but you can increase or decrease this value depending on your hardware. If a sequence's length exceeds the threshold, the script will split the sequence into smaller segments with overlapping segments (length 200) and embed each sub-sequence separately. The overlaps are averaged to produce the final embedding. This method of embedding sequences gave similar results to embedding the full sequence, but is significantly faster and allows for the embedding of any sequence regardless of it's length.

If you have access to one or more GPUs, you can also use the `--gpu` option to specify the number of GPUs to embed sequences with, otherwise it will default to using the CPU. Generating fingerprints from these embeddings is much less memory intensive, but is performed only on the CPU. You can specify the number of CPU cores to use with the `--cpu` parameter to fingerprint multiple sequences at once (default is 1). For example, the command below will embed sequences on one GPU and fingerprint them on 12 CPU cores:

```
python src/make_db.py --fafile <.fa file> --dbfile <output> --gpu 1 --cpu 12
```

## Comparing DCT fingerprints
Once the DCT fingerprints are generated, they can be used for database searching by running the query.db script. As an example, we will use the files under test/ for demonstration purposes. First run the following commands to generate the fingerprint database:

similarity detection by running the dct-sim.py script. This script can be run in three different modes. We will use the examples under test/ for demonstration purposes.

```
cd test
python ../src/make_db.py --fafile example.fasta --dbfile example
```

### Database search
To search a query fasta file/database file against a target database, run the following command:

```
python ../src/query_db.py --query example.fasta --db example.db --out example-search.txt --khits 50
```

This will generate a file containing the top 50 hits between each query sequence and the target database. Each hit lists the query sequence and its predicted domain region, the target sequence and predicted domain region, and the similarity score (1 minus the normalized L1 distance) between the two fingerprints. The scores range between 0 and 1, with 0 indicating no similarity and 1 highest similarity. Typically, DCTdomain score of 0.25 or DCTglobal score of 0.20 indicates the two proteins share some similarity. 


Like in src/make_db.py, you can change the `--maxlen`, `--cpu`, and `--gpu` parameters to suit your system's hardware (if querying a fasta file).

### Pairwise comparision for a given list of protein pairs
If you have a list of specific protein pairs that you want to compare, you can use the dct-sim.py script. The file of pairs should contain two columns, with the first column containing the ID of the first protein and the second column containing the ID of the second protein. Check example.pair for an example. The third 'label' column is not necessary to run the script.

```
python ../src/dct-sim.py --pair example.pair --dct example-dct.npz --output example-dctsim.txt 
```

The output file contains rows listing the similarity between a pair of proteins, with the first two columns showing the IDs of the proteins, followed by two similarty scores of the DCT fingerprints DCTdomain and DCTglobal.

Also see examples under bench/homo benchmark. 

## Benchmarks
See benchmarks and corresponding readmes under bench/ folder.

domcut -- domain segmentation on FUpred benchmark

G7PD -- pairs of G6PD containing proteins

homo -- homology detection benchmarks

cathdb -- database searching on CATH20

## Final words about RecCut:

If RecCut doesn't work properly on your computer, you may consider to re-compile it as follows:

```
cd src/
g++ -o RecCut RecCut.cpp
```
