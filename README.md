# What's here

- `run.py` - script to run bespokefit on input molecule. Takes a SMILES string and ID. Creates a folder with results based on the input ID. Command: `python run.py -smiles 'BrCO' -id BrCO`
- `run_all.py` - Script to run `run.py` on a CSV containing multiple SMILES. Command: `python run_all.py -fragments_csv_file host-systems/fragments.csv`
- `How_to_train_your_forcefield.ipynb` - Notebook taken from OpenFF blogpost on how to use Bespokefit.
- `host-systems/` - Contains original host files and, fragments of the host, images of the fragments, and a `.csv` file that contains fragment ID's and SMILES strings.
- `OA-fragment-*/` - Bespokefit results for a specific host fragment. 
