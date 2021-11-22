import argparse
import os

import numpy as np
from tqdm import tqdm

# predictions indices
from utils import parse_preds_file

INDEX_INDEX = 0
METAL_PROB = 1
METAL_CLASS = 2
NUCLEAR_PROB = 3
NUCLEAR_CLASS = 4
SMALL_PROB = 5
SMALL_CLASS = 6
ANY_CLASS = 7

def adjust_prob(probs : list(), distances, valid_distances, C : int) -> list():

    if not valid_distances.any():
        return probs

    boni = []
    for index in range(len(probs)):

        # get distances from index
        distance_from_index = distances[index]

        # get valid layer
        valid_from_index = valid_distances[index].astype(int)

        # ignore position to itself by setting valid layer to 0
        valid_from_index[index] = 0

        valid_distances_from_index = distance_from_index[valid_from_index.nonzero()]
        if len(valid_distances_from_index) == 0:
            boni.append(np.nan)
            continue
        valid_probs = probs[valid_from_index.nonzero()]
        boni.append(C * np.mean(np.array(valid_probs) / np.array(valid_distances_from_index)))

    boni = np.array(boni)

    # if index didnt have any valid distances, we set bonus to Nan. Now we change this to the closest recorded bonus from index
    no_bonus_indices = np.where(np.isnan(boni))[0]
    with_bonus = np.where(np.isnan(boni) == False)[0]
    for index in no_bonus_indices:
        index_2 = np.argmin(np.abs(with_bonus - index))
        boni[index] = boni[with_bonus[index_2]]

    return np.array(probs) + boni


def adjust_probs(small_probs, metal_probs, nuclear_probs, distances, valid_distances=None, C : int = 20):
    if valid_distances is None:
        valid_distances = np.ones(shape=(len(distances[0]), len(distances[0])))

    small_probs = np.array(small_probs)
    metal_probs = np.array(metal_probs)
    nuclear_probs = np.array(nuclear_probs)

    new_small_probs = adjust_prob(small_probs, distances, valid_distances, C=C)
    new_metal_probs = adjust_prob(metal_probs, distances, valid_distances, C=C)
    new_nuclear_probs = adjust_prob(nuclear_probs, distances, valid_distances, C=C)

    return new_small_probs, new_metal_probs, new_nuclear_probs


# writes adjusted probs in the same formet bindPredict does
def write_to_file(filename, outdir, modified_small_probs, modified_metal_probs, modified_nuclear_probs, cutoff = 0.5):
    with open(os.path.join(outdir, filename), "w") as file:
        file.write("Position\tMetal.Proba\tMetal.Class\tNuclear.Proba\tNuclear.Class\tSmall.Proba\tSmall.Class\tAny.Class\n")
        for index in range(len(modified_metal_probs)):
            file.write(f"{index + 1}\t{modified_metal_probs[index]:.3f}\t{'b' if modified_metal_probs[index] > cutoff else 'nb'}"
                       f"\t{modified_nuclear_probs[index]:.3f}\t{'b' if modified_nuclear_probs[index] > cutoff else 'nb'}"
                       f"\t{modified_small_probs[index]:.3f}\t{'b' if modified_small_probs[index] > cutoff else 'nb'}"
                       f"\t{'b' if modified_small_probs[index] > cutoff or modified_metal_probs[index] > cutoff or modified_nuclear_probs[index] > cutoff else 'nb'}\n")


def main(outdir : str, preds_dir : str, distance_maps_dir : str, C : int = 20, layer_index: int = 0):
    for filename in tqdm(os.listdir(preds_dir)):
        try:
            uniprot_id = filename.split(".")[0]
            if not os.path.isfile(f"{distance_maps_dir}/{uniprot_id}.npy"):
                continue

            small_probs,  metal_probs, nuclear_probs = parse_preds_file(os.path.join(preds_dir, filename))
        except:
            print(f'error parsing binding probabilities of file: {filename}')
            continue

        distance_map = np.load(f"{args.distancemap}/{uniprot_id}.npy")
        distances = distance_map[:, :, layer_index]
        valid_layer = distance_map[:, :, 4]

        modified_small_probs, modified_metal_probs, modified_nuclear_probs = adjust_probs(small_probs=small_probs, metal_probs=metal_probs, nuclear_probs=nuclear_probs, distances=distances, valid_distances=valid_layer, C=C)

        write_to_file(filename=f'{uniprot_id}.bindAdjust_out', outdir=outdir, modified_small_probs=modified_small_probs, modified_metal_probs=modified_metal_probs, modified_nuclear_probs=modified_nuclear_probs)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='bindAdjust modifies the binding probabilities of protein residues by taking into account the probability of each residue and the distances between them. The tool requires protein binding probabilites and distance maps to function.')

    parser.add_argument('-p', '--predsdir', required=True, help='directory containing predictions in specific format, see sample file. If your predictions are not available in this specific format. Please use the functions directly.')
    parser.add_argument('-o', '--outdir', required=True, help='output directory')
    parser.add_argument('-d', '--distancemap', required=True, help='directory containing protein distance maps, see sample file of distance map for required file structure and name.')

    parser.add_argument('-c', '--coef', required=False, type=int, default=20, help='coefficient used in bindAdjust. The larger C, the larger the modification.')
    parser.add_argument('-l', '--layer', required=False, type=int, default=3, help='index of distance map layer. Options: 0 N, 1 C-alpha, 2 C-beta and 3 backbone C distances')

    args = parser.parse_args()


    if not os.path.isdir(args.outdir):
        print('provided output directory does not exist or is not a directory.')
        exit(0)

    if not os.path.exists(args.distancemap):
        print('provided distance map directory does not exist or is not a directory.')
        exit(0)

    if not os.path.isdir(args.predsdir):
        print('provided prediction directory does not exist or is not a directory.')
        exit(0)

    main(outdir = args.outdir, preds_dir = args.predsdir, distance_maps_dir = args.distancemap, C = args.coef, layer_index=args.layer)