# gnina_tutorial
# GNINA docking tutorial

This tutorial shows common GNINA docking workflows: redocking, docking from a random conformer, blind docking, flexible docking, scoring/rescoring, and a basic virtual screening example. Commands below follow the GNINA software (Roes lab) — see the paper: https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00522-2

Prerequisites
- GNINA installed and available in your PATH (or use full binary path, e.g. /usr/local/gnina/gnina.1.3.2).
- OpenBabel (obabel, obrms) installed.
- wget, gzip, awk, python3 with matplotlib and pandas for plotting.
- A molecular viewer (VMD, PyMOL, ChimeraX) to inspect results.

Repository layout used in this tutorial (created on the fly):
- re-docking/         — reproduce bound pose (redocking)
- docking/            — docking from random conformer
- rescoring/          — scoring-only experiments
- flexible/           — flexible docking experiments
- VS/                 — virtual screening example
- run_dock_exhaustiveness.sh — helper to test exhaustiveness
- plot_exhaustiveness.py     — plot RMSD vs exhaustiveness

Quick notes about common commands used here
- gnina -r <receptor> -l <ligand> [options] -o <out>
- --autobox_ligand <lig> : center the docking box on a ligand (e.g. the bound ligand)
- --exhaustiveness N : increase search thoroughness (default: 8)
- --score_only : run scoring only (no pose generation)
- --flexres A:52,A:103 : specify flexible residues
- --flexdist N and --flexdist_ligand <lig> : set distance-based flexible region
- obrms -firstonly ref.pdb poses.pdb : compute RMSD of the first (top) pose vs reference

1) Redocking: recover the exact bound conformation
-------------------------------------------------
Purpose: show GNINA can re-create the experimental bound pose.

Commands (inside a working directory, e.g. re-docking):
1. Make a directory and download the PDB
```
   mkdir re-docking
   cd re-docking
   wget https://files.rcsb.org/download/3ERK.pdb
```
3. Extract protein atoms only (ATOM records)
   ```
   grep "^ATOM" 3ERK.pdb > rec.pdb
```
4. Fix connectivity (OpenBabel will add bond info)
```
   obabel rec.pdb -O rec-fix.pdb
```
6. Extract ligand (SB4 is the ligand code in 3ERK)
```
   grep SB4 3ERK.pdb > lig.pdb
```
7. Perform a simple docking (autobox_ligand centers the box on the bound ligand)
```
   gnina -r rec.pdb -l lig.pdb --autobox_ligand lig.pdb --seed 0 -o docked.pdb
```
8. Compute RMSD of top pose vs experimental ligand
   ```
   obrms -firstonly lig.pdb docked.pdb
```
Open docked.pdb in VMD / PyMOL and compare poses visually.

Assignment 1 (suggested)
- Investigate the effect of --exhaustiveness on the docking RMSD.
- Use the helper script run_dock_exhaustiveness.sh included in this repo to automate runs and produce a CSV of RMSD vs exhaustiveness values.
- Plot RMSD vs exhaustiveness with plot_exhaustiveness.py and determine the lowest exhaustiveness that yields converged (stable/acceptable) RMSD. Default is 8.

2) Docking from a random conformer (same system, new starting ligand)
--------------------------------------------------------------------
Purpose: check robustness to ligand starting conformation.

1. Go up one directory and start a new folder
   ```
   cd ..
   mkdir docking
   cd docking
```
2. Copy receptor and reference ligand for autoboxing
   ```
cp ../re-docking/rec.pdb ../re-docking/lig.pdb .
```
3. Generate a random 3D conformer from SMILES using OpenBabel
   ```
obabel -:'C1CNCCC1N2C=NC(=C2C3=NC(=NC=C3)N)C4=CC=C(C=C4)F' -O lig-random.sdf --gen3D
```
4. Dock the random conformer (example uses exhaustiveness 8)
   ```
gnina -r ../rec.pdb -l lig-random.sdf --autobox_ligand ../re-docking/lig.pdb --seed 0 -o docked_random.pdb --exhaustiveness 8
```
5. RMSD of the best pose vs experimental ligand
   ```
obrms -firstonly ../re-docking/lig.pdb docked_random.pdb
```
Investigate: which starting conformation gives better RMSD and CNN score (GNINA prints CNN scoring information). Compare affinity (CNN score, Vina/Vinardo scores) as reported in the output SDF/PDB.

3) Scoring-only / rescoring
--------------------------
Purpose: evaluate different scoring functions on fixed poses.

- Score only with default scoring:
  ```
gnina --score_only -r ../re-docking/rec.pdb -l ../re-docking/lig.pdb --verbosity=2
```
- List available scoring functions:
  ```
gnina --help | grep scoring | head -3
```
- Score with classical scoring (example: ad4_scoring):
  ```
gnina --score_only -r rec.pdb -l lig.pdb --verbosity=2 --scoring ad4_scoring
```
- Print scoring terms:
  ```
gnina --score_only -r rec.pdb -l lig.pdb --verbosity=2 --scoring ad4_scoring --print_terms
```
- Custom scoring terms file:
  ```
gnina --score_only -r rec.pdb -l lig.pdb --score_only --custom_scoring scoring-terms-all.txt
```
Try the convolutional neural network (CNN) scoring family:
- See available CNN models:
  gnina --help | grep "cnn arg" -A 12

- Score with CNN (CNN score is between 0 and 1; higher is better)
  ```
gnina --score_only -r rec.pdb -l lig.pdb | grep CNN
```
4) Blind docking
----------------
Purpose: search the whole receptor surface for possible binding sites.

- Blind docking by centering autobox on receptor (autobox_ligand rec.pdb)
  ```
gnina -r rec.pdb -l lig.pdb --autobox_ligand rec.pdb -o docked-blind.pdb --seed 0
```
- For a random ligand conformer:
  ```
obrms -firstonly lig.pdb docked-blind.pdb | head -n 10
```
Notes:
- Blind docking is slower and may require larger exhaustiveness.
- Evaluate top clusters/poses and their CNN scores to prioritize predicted sites.

5) Flexible docking
-------------------
Purpose: allow receptor residues (sidechains / local backbone) to move during docking.

1. Dock ligand from 3ERK into a similar protein (4ERK) without flexibility:
   ```
wget http://files.rcsb.org/download/4ERK.pdb
```
```
grep "^ATOM" 4ERK.pdb > rec2.pdb
obabel rec2.pdb -O rec2-fix.pdb
grep OLO 4ERK.pdb > lig2.pdb
```
```
gnina -r rec2-fix.pdb -l ../re-docking/lig.pdb --autobox_ligand lig2.pdb --seed 0 -o docked-lig-onto-4ERK.pdb
```

```
obrms -firstonly ../re-docking/lig.pdb docked-lig-onto-4ERK.pdb
```
2. Try higher exhaustiveness if initial docking is poor:
   ```
gnina -r rec2-fix.pdb -l ../re-docking/lig.pdb --autobox_ligand lig2.pdb --seed 0 -o docked-lig-onto-4ERK.pdb --exhaustiveness 512
```
3. Enable flexible residues by distance from the autobox ligand:
   ```
gnina -r rec2-fix.pdb -l ../re-docking/lig.pdb --autobox_ligand lig2.pdb --seed 0 -o docked-lig-onto-4ERK-flex.pdb --flexdist 4 --flexdist_ligand lig2.pdb --out_flex flexout.pdb
```
4. Or specify exact residues to be flexible:
   ```
gnina -r rec2-fix.pdb -l ../re-docking/lig.pdb --autobox_ligand lig2.pdb --seed 0 -o docked-lig-onto-4ERK-flex2.pdb --exhaustiveness 16 --flexres A:52,A:103 --out_flex flexout2.pdb
```
Assignment 2 (suggested)
- Investigate convergence of results vs. exhaustiveness in the flexible and rigid cases and compare RMSD / CNN score.
- Which flexible setup (distance-based vs explicit residues) gives the best recovery for the 4ERK system?

6) Virtual screening (VS) example
--------------------------------
Purpose: score a ligand database against a receptor and assess enrichment.

1. Example files used below:
```
wget http://files.rcsb.org/download/4PPS.pdb
```

```
wget http://bits.csb.pitt.edu/files/workshop_minimized_results.sdf.gz
```

```
   grep ^ATOM 4PPS.pdb > errec.pdb
```
2. A simple VS run with VINARDO scoring:
   ```
gnina -r errec.pdb -l workshop_minimized_results.sdf.gz --minimize -o gnina_scored_vinardo.sdf.gz --scoring vinardo
```
3. Compute ROC / AUC:
   - The file how-to-run-VS1.txt contains instructions and pointers to VS1.py (script included in this repo or workshop). Use that script to get ROC/AUC and typical enrichment metrics.

Assignment 3 (suggested)
- Try other scoring functions (cnn family, classic scorings) and compare AUC/enrichment.
- Try rescoring GNINA poses with CNN to see improvement vs classical scoring.

Helper scripts included
- run_dock_exhaustiveness.sh : run GNINA for a list of exhaustiveness values and collect RMSD of the top pose into a CSV.
- plot_exhaustiveness.py : read CSV and produce a plot (matplotlib).

Final notes and tips
- Use --seed for reproducibility.
- Increase --exhaustiveness for more thorough search; the runtime increases roughly linearly with it.
- Use --cnn_only if you only want CNN-based pose ranking (check your GNINA build options).
- When comparing RMSDs, make sure you align atoms consistently (obrms -firstonly is useful) and be aware of symmetric ligands.

References
- GNINA paper: Roes et al., Journal of Cheminformatics (2021). https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00522-2

If you want, I can:
- Add the helper scripts to the repository (run_dock_exhaustiveness.sh and plot_exhaustiveness.py), or
- Generate a ready-to-run example that executes all steps for 3ERK (requires GNINA/OpenBabel installed on your machine).
