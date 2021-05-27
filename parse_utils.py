import os
import time
import imageio as imageio
from pymol import cmd
import argparse

# ligand map indices
PDB_ID_INDEX = 0
RES_INDEX = 2
UNIPROT_ID_INDEX = 3

# predictions indices
INDEX_INDEX = 0
METAL_PROB = 1
METAL_CLASS = 2
NUCLEAR_PROB = 3
NUCLEAR_CLASS = 4
SMALL_PROB = 5
SMALL_CLASS = 6
ANY_CLASS = 7

def parse_uniprot_ids(filepath: str) -> set():
    ids = list()
    with open(filepath, 'r', ) as file:
        for line in file:
            ids.append(line.strip())
    return ids


def parse_ligand_map(filepath: str):
    uniprot_to_best_res = {}
    uniprot_to_offset = {}
    uniprot_to_pdb ={}
    with open(filepath, "r") as map_file:
        map_file.readline()  # skip header line
        for line in map_file:
            stripped_line = line  # strip white spaces
            values = stripped_line.split("\t")  # split at tab
            uniprot_id = values[UNIPROT_ID_INDEX]
            resolution = float(values[RES_INDEX])

            # check if we have this id already or if we have higher resolution structures saved
            if uniprot_id not in uniprot_to_best_res.keys() or resolution < uniprot_to_best_res[uniprot_id]:
                uniprot_to_pdb[uniprot_id] = values[PDB_ID_INDEX]
                uniprot_to_best_res[uniprot_id] = resolution

            if uniprot_id not in uniprot_to_offset.keys():
                uniprot_to_offset[uniprot_id] = int(values[-2].split(",")[0]) - int(values[-1].split(",")[0])
    return uniprot_to_pdb, uniprot_to_offset


def parse_trues(filepath: str) -> dict():
    # iterate over true annotations metal for now
    uniport_to_true_annots = {}
    with open(filepath, "r") as annot_file:
        for line in annot_file:
            stripped_line = line.strip()  # strip white spaces
            values = stripped_line.split("\t")  # split at tab
            uniport_id = values[0]
            uniport_to_true_annots[uniport_id] = set(values[1].split(","))

    return uniport_to_true_annots


def parse_preds(filedir: str, threshold: float):
    #TODO: add parameter for ligand type, nulcear hardcoded rn
    uniprot_to_preds_dict = {}
    # iterate over files
    for filename in os.listdir(filedir):
        # print(filename)
        uniprot_id = filename.split(".")[0]
        with open(os.path.join(filedir, filename), "r") as predictions:
            predictions.readline()  # skip first line
            for line in predictions:
                values = line.split()
                if float(values[SMALL_PROB]) > threshold:
                    if uniprot_id in uniprot_to_preds_dict.keys():
                        uniprot_to_preds_dict[uniprot_id].add(values[INDEX_INDEX])
                    else:
                        uniprot_to_preds_dict[uniprot_id] = set()
                        uniprot_to_preds_dict[uniprot_id].add(values[INDEX_INDEX])

    return uniprot_to_preds_dict