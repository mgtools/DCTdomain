# DCTdomain
Protein domain-level embedding via ESM-2 and DCT for homology detection and similarity search.

Paper: https://www.biorxiv.org/content/10.1101/2023.11.27.567555v2

## Approach
An overview of the pipeline is visualized in the figure below:

![This is an image](https://github.com/mgtools/DCTdomain/blob/main/misc/DCTdomain-diag.png)

## Programs/scripts
This package contains programs and scripts for the following tasks:
1) Programs/scripts required for generating ESM-2 embeddings and contact maps
2) Apply RecCut to predict the domains given contact maps
3) Generate DCT fingerprints (DCTdomain and DCTglobal) given predicted domains
4) Calculate similarity between DCT fingerprints 

## Generating DCT Fingerprints
The pipeline to generate DCT fingerprints can be called by running the bash script as shown below:

```
bash src/fingerprint.sh <.fa file> <output1-name>
```

fingerprint.sh takes a fasta file with any number of sequences. This will produce three output files, output.list which contains a list of all embedded protein sequence, output.dom which contains a list of predicted domains for each protein in the list file, and output-dct.npz which contains one (single-domain proteins) or more DCT fingerprints (multi-domain proteins) for each protein. Each individual embedding and top contacts for each sequence are saved under separate directories in case the embedding process is interrupted.

### Embedding protein sequences
Depending on the length of the protein sequences in your fasta file, you may need to adjust the '--maxlen' parameter when calling embed-mult.py in fingerprint.sh. The default value is 1000, for which a CPU with 32GB of RAM should be sufficient. You can increase or decrease this value depending on your hardware. If you have a GPU, you can also use the '--gpu' option to speed up the embedding process.

If a sequence's length exceeds the threshold, the script will split the sequence into smaller segments with overlapping segments (length 200) and embed each sub-sequence separately. The overlaps are averaged to produce the final embedding. This method of embedding sequences gave similar results to embedding the full sequence, but is significantly faster and allows for the embedding of any sequence regardless of it's length.

There is also a '--checkpoint' parameter that can be used to specify which ESM-2 model to embed with. The default is esm2_t30_150M_UR50D, which is the model used for the results found in the paper. Contact maps from t33 were found be more accurate for domain prediction (results in Supplementary Information), but using t30 was found to be more accurate for our homology inference task. You can use esm2_t33_650_UR50D by specifying:

```
python src/embed-mult.py --input <input.fa> --npy <directory> --checkpoint t33
```

## Comparing DCT fingerprints
Once the DCT fingerprints are generated, they can be used for similarity detection by running the dct-sim.py script. This script can be run in three different modes. We will use the examples under test/ for demonstration purposes.

```
cd test
```

### All against all pairwise comparison
By default, dct-sim.py will perform an all-against-all comparison of the proteins in the .npz file. You can capture the output by using the --output option.

```
python ../src/dct-sim.py --dct example-dct.npz --output example-all-vs-all-dctsim.txt
```

### Pairwise comparision for a given list of protein pairs
If you have a list of specific protein pairs that you want to compare, you can use the --pair option. The file of pairs should contain two columns, with the first column containing the ID of the first protein and the second column containing the ID of the second protein. Check example.pair for an example. The third 'label' column is not necessary to run the script.

```
python ../src/dct-sim.py --pair example.pair --dct example-dct.npz --output example-dctsim.txt 
```

Also see examples under bench/homo benchmark. 

### All against database
In this example we use the same npz as the query and as the database. All proteins given in --dct will be compared against all proteins given in --db. 

```
python ../src/dct-sim.py --dct example-dct.npz --db example-dct.npz --output example-all-vs-db-dctsim.txt
```

### Results and DCT fingerprint similarty scores

In the main output file, each row shows the similarity between a pair of proteins, with the first two columns showing the IDs of the proteins, followed by two similarty scores of the DCT fingerprints DCTdomain and DCTglobal. The scores range between 0 and 1, with 0 indicating no similarity and 1 highest similarity. Typically, DCTdomain score of 0.25 or DCTglobal score of 0.20 indicates the two proteins share some similarity. 

## Benchmarks
See benchmarks and corresponding readmes under bench/ folder.

domcut -- domain segmentation on FUpred benchmark

G7PD -- pairs of G6PD containing proteins

homo -- homology detection benchmarks

## Requirements
To install DCTdomain, clone this repository:

```
git clone https://github.com/mgtools/DCTdomain
```

And install the following dependences:

python3 >= 3.9
python packages:

    -fair-esm
    -torch
    -numpy
    -scipy

Or you can create and activate a conda environment as follows:

```
cd DCTdomain/
conda create -f env.yml
conda activate dctdomain
```

## Final words about RecCut:

If RecCut doesn't work properly on your computer, you may consider to re-compile it as follows:

```
cd src/
g++ -o RecCut RecCut.cpp
```
