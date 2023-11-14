# DCTdomain
Protein domain-level embedding via ESM-2 and DCT for homology detection and similarity search.

## Approach
An overview of the pipeline is visualized in the figure below:

![This is an image](https://github.com/mgtools/DCTdomain/blob/main/misc/DCTdomain-diag.png)

## Programs/scripts
this package contains programs and scripts for the different tasks
1) all programs/script required for generating ESM-2 embedding and contact maps
2) apply RecCut to disect the domains given contact map
3) generate DCT fingerprints (DCTdomain and DCTglobal) given predicted domains
4) applications of DCT fingerprints 

## Running Pipeline
This pipeline is implemented in two scripts, a bash script that generates the DCT fingerprints and a python script that calculates the similarity between them. You can run them as shown below:

```
bash src/fingerprint.sh <.fa file> <.list file> <output1-name>
python src/dct-sim.py --pair <.pair file> --dct <output1-name>-dct.npz --output <output2-name>
```

fingerprint.sh takes a fasta file with any number of sequences and file with a list of each sequence you would like a fingerprint for (see files under test/ for examples). This will produce two output files, <output>.dom which contains a list of predicted domains for each protein in the fasta file, and <output>-dct.npz which contains a DCT fingerprint for each predicted domain.

dct-sim.py takes a file with a list of each sequence pair you would like to calculate the distance between and a .npz file containing DCT fingerprints from src/fingerprint.sh (again see files under test/ for examples). This will produce a file with the distances between each pair of sequences.

## Using DCT fingerprints in three different modes
The script dct-sim.py can be run in three different modes. 

We will sse the examples under test for demonstration purposes.

cd test

### Pairwise comparision for a given list of protein pairs (--pair)
python ../src/dct-sim.py --pair example.pair --dct example-dct.npz --output example-dctsim.txt 

### All against all pariwise comparison for all the proteins
python ../src/dct-sim.py --dct example-dct.npz --output example-all-vs-all-dctsim.txt 

In this example, all proteins with DCT fingerprints in example-dct.npz will be compared against each other.

### All against database
python ../src/dct-sim.py --dct example-dct.npz --db example-dct.npz --output example-all-vs-db-dctsim.txt 

In this example we use the same npz as the query and as the database. All proteins given in --dct will be compared against all proteins given in --db. 

## Benchmarks
See benchmarks and corresponding readmes under bench/ folder.

domcut -- domain segmentation on FUpred benchmark

G7PD -- pairs of G6PD containing proteins

homo -- homology detection benchmarks

## Requirements
To install DCTDomain, clone this repository:

```
git clone https://github.com/mgtools/DCTdomain
```

And install the following dependences:

python3
python packages:
    -fair-esm
    -torch
    -numpy
    -scipy
