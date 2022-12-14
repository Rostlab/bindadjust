import os

# predictions indices
INDEX_INDEX = 0
METAL_PROB = 1
METAL_CLASS = 2
NUCLEAR_PROB = 3
NUCLEAR_CLASS = 4
SMALL_PROB = 5
SMALL_CLASS = 6
ANY_CLASS = 7

# ligand map indices
PDB_ID_INDEX = 0
RES_INDEX = 2
UNIPROT_ID_INDEX = 3


def parse_trues(filepath: str) -> dict():
    try:
        uniport_to_true = {}
        with open(filepath, "r") as file:
            for line in file:
                values = line.strip().split("\t")
                uniport_id = values[0]
                uniport_to_true[uniport_id] = set([int(value) - 1 for value in values[1].split(",")])
        return uniport_to_true
    except:
        print('error parsing true file, pls check if file has correct format')

def parse_preds_file(filepath : str):
    small_probs = []
    metal_probs = []
    nuclear_probs = []

    with open(filepath, "r") as predictions:
        predictions.readline()  # skip header line
        for line in predictions:
            values = line.split()
            small_probs.append(float(values[SMALL_PROB]))
            metal_probs.append(float(values[METAL_PROB]))
            nuclear_probs.append(float(values[NUCLEAR_PROB]))

    return small_probs, metal_probs, nuclear_probs



def parse_preds(filedir: str, ligand_type: str) -> dict():
    uniprot_to_preds_dict = {}
    for filename in os.listdir(filedir):
        try:
            uniprot_id = filename.split(".")[0]
            small_probs, metal_probs, nuclear_probs = parse_preds_file(os.path.join(filedir, filename))
            if ligand_type == 'small': uniprot_to_preds_dict[uniprot_id] = small_probs
            elif ligand_type == 'metal': uniprot_to_preds_dict[uniprot_id] = metal_probs
            elif ligand_type == 'nuclear': uniprot_to_preds_dict[uniprot_id] = nuclear_probs
        except:
            print(f'error parsing prediction {filename = }. Pls make sure no other files are in predictions folder.')

    return uniprot_to_preds_dict


def parse_ligand_map(filepath: str):
    try:
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

                if uniprot_id not in uniprot_to_best_res.keys() or resolution < uniprot_to_best_res[uniprot_id]:
                    uniprot_to_pdb[uniprot_id] = values[PDB_ID_INDEX]
                    uniprot_to_best_res[uniprot_id] = resolution
                    uniprot_to_offset[uniprot_id] = int(values[-2].split(",")[0]) - int(values[-1].split(",")[0])
    except:
        print('error parsing ligand map, pls use correct format or use stored ligand map')

    return uniprot_to_pdb, uniprot_to_offset