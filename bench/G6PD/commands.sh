#Inputs: G6PD.fasta, G6PD.info, G6PD.pair

#Compute DCTdomain fingerprints
bash ../../src/fingerprint.sh G6PD.fasta G6PD

#Compute similarity scores using DCTdomain fingerprints (save the results to G6PD-dctsim.txt
python ../../src/dct-sim.py --pair G6PD.pair --dct G6PD-dct.npz --output G6PD-dctsim.txt

#to prepare plots showing the comparison of the score distribions (as those in the paper)
#python hist.py
