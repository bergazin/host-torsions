# What's here
Scripts for running Bespokefit on host systems and results.

**Host to do list:**
- [ ] GDCC hosts from SAMPL7 and SAMPL8 https://github.com/samplchallenges/SAMPL7/tree/master/host_guest/GDCC_and_guests (or https://github.com/samplchallenges/sampl8) **[WIP]**
  - [ ] OA **[WIP]**
  - [ ] exo OA **[TO DO]**
- [ ] Cyclodextrins. set of modified cyclodextrin derivatives: https://github.com/samplchallenges/SAMPL7/tree/master/host_guest/cyclodextrin_derivatives **[TO DO]**
- [ ] MAYBE the SAMPL9 WP6 host: https://github.com/samplchallenges/sampl9

Setting up Bespokefit for this project:
1. Clone openFF bespokefit
2. Add `xtb-python` to `bespoke-fit/devtools/conda-envs/test-env.yaml`
2. Create environemnt from testing environment `conda env update -n bespoke --file bespoke-fit/devtools/conda-envs/test-env.yaml`. Install additonal packages: `conda install -c conda-forge mamba seaborn`.
4. Execute `python setup.py develop` in `bespoke-fit/`
5. Run `pip install git+git://github.com/openforcefield/bespoke-fit.git@931bb66bb3700208616e41bfde1a80d90a37d74f` to get desired environment.


## 📂 Directory
- [`run.py`](run.py) - script to run bespokefit on input molecule. Takes a SMILES string and ID. Creates a folder with results based on the input ID. Command: `python run.py -smiles 'BrCO' -id BrCO`
- [`run_all.py`](run_all.py) - Script to run `run.py` on a CSV containing multiple SMILES. Command: `python run_all.py -fragments_csv_file host-systems/fragments.csv`
- [`How_to_train_your_forcefield.ipynb`](How_to_train_your_forcefield.ipynb) - Notebook taken from OpenFF blogpost on how to use Bespokefit.
- [`host-systems/`](host-systems/) - Contains original host files and, fragments of the host, images of the fragments, and a `.csv` file that contains fragment ID's and SMILES strings.
- `OA-fragment-*/` - Bespokefit results for a specific host fragment.
