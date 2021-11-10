#!/usr/bin/env python3
import pandas as pd
import os
import logging
import argparse
import sys
from tqdm import tqdm


logging.getLogger().setLevel(logging.INFO)

def parse_options() -> argparse.Namespace:
    """Parse main options."""
    parser = argparse.ArgumentParser(description="Command line arguments")
    parser.add_argument('-fragments_csv_file', type=str, metavar='str', help='CSV file that contains fragments IDs and SMILES strings')
    return parser.parse_args()


def main() -> int:
    """
    Runs each molecule in input CSV through Bespokefit workflow.
    """

    args = parse_options()

    df = pd.read_csv("host-systems/fragments.csv")
    for idx in tqdm(range(len(df)), desc="Progress"):
    #for idx in range(len(df)):
        fragment_id = df.loc[idx, 'fragment_id']
        smiles = df.loc[idx, 'smiles_from_chemdraw']
        logging.warning(f"Working on: {fragment_id}")
        os.system(f"python run.py -smiles '{smiles}' -id {fragment_id}")
        logging.warning(f"Bespoke executed for: {fragment_id}")


    return 0

if __name__ == "__main__":
    sys.exit(main())
