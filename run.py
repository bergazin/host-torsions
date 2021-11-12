#!/usr/bin/env python3
from openff.bespokefit.schema.optimizers import ForceBalanceSchema
from openff.bespokefit.workflows import BespokeWorkflowFactory
from openff.fragmenter.fragment import WBOFragmenter
from openff.bespokefit.schema.targets import TorsionProfileTargetSchema
from openff.qcsubmit.common_structures import QCSpec
from openff.bespokefit.schema.smirnoff import  ProperTorsionHyperparameters
from openff.toolkit.topology import Molecule
from rdkit.Chem import Draw, rdDepictor
from openff.bespokefit.executor import BespokeExecutor, wait_until_complete
import os
import matplotlib.pyplot as plt
import seaborn as sns
#rom IPython.display import IFrame
import pandas as pd
from simtk import unit
import shutil
import logging
import argparse
import sys
import json

logging.getLogger().setLevel(logging.INFO)

def parse_options() -> argparse.Namespace:
    """Parse main options."""
    parser = argparse.ArgumentParser(description="Command line arguments")
    parser.add_argument('-smiles', type=str, metavar='str', help='SMILES string')
    parser.add_argument('-id', type=str, metavar='str', help='Name of molecule (will be used to generate output SDF and for output files)')
    return parser.parse_args()


def main() -> int:
    """
    Run bespokefit on input molecules.
    """

    args = parse_options()

    output_directory = f"./{args.id}"
    if not os.path.exists(os.getcwd() + '/' + output_directory):
        os.makedirs(os.getcwd() + '/' + output_directory, exist_ok=True)

    fb = ForceBalanceSchema()
    fb.penalty_type = "L1"

    workflow = BespokeWorkflowFactory(
        fragmentation_engine=None,
        optimizer=fb,
        parameter_hyperparameters=[],
        target_templates=[],
        default_qc_specs=[]
    )
    workflow.dict()


    fragmenter = WBOFragmenter()
    workflow.fragmentation_engine = fragmenter
    workflow.target_torsion_smirks = ["[*]~[!$(*#*)&!D1:1]-,=;!@[!$(*#*)&!D1:2]~[*]"]

    target = TorsionProfileTargetSchema()
    workflow.target_templates = [target, ]
    target.dict()


    xtb_spec = QCSpec(
        method="gfn2xtb",
        basis=None,
        program="xtb",
        spec_name="xtb",
        spec_description="gfn2xtb"
    )

    workflow.default_qc_specs = [xtb_spec]

    prior = ProperTorsionHyperparameters()
    workflow.parameter_hyperparameters = [prior]

    workflow.generate_bespoke_terms = True
    workflow.expand_torsion_terms = True

    logging.warning(f'Writing out workflow to {args.id}_workflow.json and {args.id}_workflow.yaml')
    # write to json
    workflow.export_factory(f"./{args.id}/{args.id}_workflow.json")
    # or yaml
    workflow.export_factory(f"./{args.id}/{args.id}_workflow.yaml")

    target_molecule = Molecule.from_smiles(args.smiles)
    target_molecule.generate_conformers(n_conformers=1)
    # write to file
    target_molecule.to_file(file_path=f"./{args.id}/{args.id}.sdf", file_format="sdf")
    # make the specific schema
    logging.warning(f'Making the specific schema for {args.id}')
    schema = workflow.optimization_schema_from_molecule(molecule=target_molecule)
    logging.warning(f'Done making the schema')

    parameters = schema.parameters
    rd_mol = target_molecule.to_rdkit()
    rdDepictor.Compute2DCoords(rd_mol)
    rd_mols = []
    atoms = []
    for parameter in parameters:
        matches = target_molecule.chemical_environment_matches(query=parameter.smirks)
        flat_matches = [atom for match in matches for atom in match]
        rd_mols.append(rd_mol)
        atoms.append(flat_matches)



    # set keep files to true so we can view the results
    os.environ["BEFLOW_KEEP_FILES"] = "True"

    # launch the executor
    logging.warning('Launching the executor')
    with BespokeExecutor(
            n_fragmenter_workers=1,
            n_qc_compute_workers=1,
            n_optimizer_workers=1) as executor:
        # grab the task id and wait for the task to finish
        task = executor.submit(input_schema=schema)
        result = wait_until_complete(optimization_id=task.id)

    print("Stages:")
    for stage in result.stages:
        print(stage.type, stage.status)


    opt_id = f"bespoke-executor/{result.results.input_schema.id}/optimize.tmp/torsion-0/iter_0002/plot_torsion.pdf"
    dst_path = f"./{args.id}/{args.id}_plot_torsion.pdf"
    shutil.copy(opt_id, dst_path)

    #IFrame(
    #    opt_id,
    #    width=900,
    #    height=600)

    #initial and new paramters
    #result.results.refit_parameter_values
    #result.results.input_schema.initial_parameter_values

    # Write old paramteres to json file (redundant)
    old_parameters = {torsion.smirks: {k: q._value for k, q in parameters.items()} for torsion, parameters in result.results.input_schema.initial_parameter_values.items()}
    with open(f"{args.id}/{args.id}_initial_parameter_values.json", "w") as f:
      json.dump(old_parameters, f)

    # Write new paramteres to json file (redundant)
    new_parameters = {torsion.smirks: {k: q._value for k, q in parameters.items()} for torsion, parameters in result.results.refit_parameter_values.items()}
    with open(f"{args.id}/{args.id}_refit_parameter_values.json", "w") as f:
      json.dump(new_parameters, f)

    parameter_data = []
    for i, (parameter, initial_values) in enumerate(result.results.input_schema.initial_parameter_values.items()):
        for final_parameter, final_values in result.results.refit_parameter_values.items():
            if parameter.smirks == final_parameter.smirks:
                for term in range(1, 5):
                    k_before = initial_values[f"k{term}"].value_in_unit(unit.kilocalorie_per_mole)
                    k_after = final_values[f"k{term}"].value_in_unit(unit.kilocalorie_per_mole)
                    parameter_data.append([parameter.smirks, f"smirks_{i}_k{term}", k_before, k_after, k_after - k_before])

    # make a pandas dataframe
    df = pd.DataFrame(parameter_data, columns=["smirks", "parameter", "before", "after", "change"])
    print("Changes:")
    print(df)
    df.to_csv(f"{args.id}/{args.id}_parameter_changes.csv")


    #%matplotlib inline

    logging.warning('Plotting parameters changes')
    plt.rc('font', size=16)
    ax = plt.figure(figsize=(5,8))
    #add start points
    ax = sns.stripplot(data=df,
                       x='before',
                       y='parameter',
                       orient='h',
                       order=df['parameter'],
                       size=10,
                       color='black')
    ax.grid(axis='y', color='0.9')
    #add arrows to plot only if the parameter changed by more than 1e-3 kcal/mol
    for i in range(len(df.index)):
        term = df.iloc[i]
        if abs(term["change"]) > 1e-3:
            if term["after"] > term["before"]:
                arrow_color = '#347768'
            elif term["after"] < term["before"]:
                arrow_color = 'red'
            else:
                arrow_color = 'black'
            ax.arrow(term["before"],
                     i,
                     term["change"],
                     0,
                     head_width=0.3,
                     head_length=0.2,
                     width=0.1,
                     fc=arrow_color,
                     ec=arrow_color)
    ax.axvline(x=0, color='0.9', ls='--', lw=2, zorder=0)
    ax.set_xlabel('Force Constant (kcal/mol)')
    ax.set_ylabel('Torsion Parameter')
    sns.despine(left=True, bottom=True)
    plt.savefig(f"{args.id}/{args.id}_parameter_change.pdf", bbox_inches = "tight")

    # copy run script to results dir too
    shutil.copy(f"./run.py", f"./{args.id}/{args.id}_run_copy.py")

    logging.warning(f'Bespoke run for {args.id} is done')
    logging.warning(f'Check the "{args.id}" folder for output files')

    logging.warning('Highlighting torsions')
    img = Draw.MolsToGridImage(rd_mols, highlightAtomLists=atoms)
    img.save(f"./{args.id}/{args.id}.png")#,bbox_inches=’tight’)
    return 0

if __name__ == "__main__":
    sys.exit(main())
