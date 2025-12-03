#!/bin/bash

for e in 8 16 24 32; do
    mkdir -p "$e"
    cd "$e"

    ln -sf ../lig.pdb .
    ln -sf ../rec.pdb .

    
    /usr/local/gnina/gnina.1.3.2 -r rec.pdb -l lig.pdb \
        --autobox_ligand lig.pdb \
        --seed 0 \
        -o docked.pdb \
        --exhaustiveness "$e"

    obrms -firstonly lig.pdb docked.pdb > rmsd.dat

    cd ..
done
