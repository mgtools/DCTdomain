#pfam-max50 bench
python ../../src/dct-sim.py --dct pfam_max50-dct.npz --pair pfam_max50.pair --output pfam_max50-dctsim.txt

#pfam-nomax50 bench
python ../../src/dct-sim.py --dct pfam_nomax50-dct.npz --pair pfam_nomax50.pair --output pfam_nomax50-dctsim.txt

#pfam-local bench
python ../../src/dct-sim.py --dct pfam_nomax50-dct.npz --pair pfam_localpfam_nomax50.pair --output pfam_localpfam_nomax50-dctsim.txt
