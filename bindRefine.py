import argparse
import numbers
import os
import warnings

import numpy as np
from tqdm import tqdm

from utils import parse_preds_file


def refine(probs, distances, valid_layer, epsilon=None, k=None, l=1, overlap=True):
    if epsilon is None and k is None:
        print("Please specify either k or epsilon. For more information on parameter usage, please read the documentation at [link] or print help using -h")
        exit(-1)

    if epsilon is not None and k is not None:
        print("Cannot use epsilon and k. Please specify only one: k or epsilon. For more information on parameter usage, please read the documentation at [link] or print help using -h")
        exit(-1)

    if epsilon is not None and (not isinstance(epsilon, numbers.Number) or epsilon < 0):
        print("epsilon must be a positive number. For more information on parameter usage, please read the documentation at [link] or print help using -h")
        exit(-1)

    if k is not None and (not type(k) == int or k < 0):
        print("k must be a positive integer. For more information on parameter usage, please read the documentation at [link] or print help using -h")
        exit(-1)

    probs = np.array(probs)
    highest_average_binding = 0
    residues_in_hotspot = None

    for center_index in range(len(probs)):
        relevant_distances = distances[center_index]

        # neighborhood definiton differs for eps and k version
        if epsilon is not None:
            # select all residues which have valid distances to center_index and distacne is below < epsilon
            valid_and_in_range = [resi_index for resi_index in range(len(relevant_distances)) if relevant_distances[resi_index] < epsilon and valid_layer[center_index][resi_index] == 1]
        elif k is not None:
            # select k-1 closest residues to center_residue
            valid_enumerated = [value for value in enumerate(relevant_distances) if valid_layer[center_index][value[0]] == 1]
            valid_enumerated.sort(key=lambda x: x[1])
            valid_and_in_range = [value[0] for value in valid_enumerated[:k]]

        # we suppress warning if valid and in range is empty
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            average_binding_prob = np.mean(probs[valid_and_in_range])

        if average_binding_prob > highest_average_binding:
            highest_average_binding = average_binding_prob
            residues_in_hotspot = valid_and_in_range

    output = np.zeros(len(probs))
    output[residues_in_hotspot] = 1
    return output


def write_to_file(filename, outdir, small_hotspot, metal_hotspot, nuclear_hotspot, small_probs, metal_probs, nuclear_probs, cutoff = 0.5):
    with open(os.path.join(outdir, filename), "w") as file:
        file.write("Position\tMetal.Proba\tMetal.Class\tMetal.Refine\tNuclear.Proba\tNuclear.Class\tNuclear.Refine\tSmall.Proba\tSmall.Class\tSmall.Refine\tAny.Class\n")
        for index in range(len(small_probs)):
            file.write(f"{index + 1}\t{metal_probs[index]:.3f}\t{'b' if metal_probs[index] > cutoff else 'nb'}\t{int(metal_hotspot[index])}"
                       f"\t{nuclear_probs[index]:.3f}\t{'b' if nuclear_probs[index] > cutoff else 'nb'}\t{int(nuclear_hotspot[index])}"
                       f"\t{small_probs[index]:.3f}\t{'b' if small_probs[index] > cutoff else 'nb'}\t{int(small_hotspot[index])}"
                       f"\t{'b' if small_probs[index] > cutoff or metal_probs[index] > cutoff or nuclear_hotspot[index] > cutoff else 'nb'}\n")


def main(outdir: str, preds_dir: str, distance_maps_dir: str, layer_index: int = 3, k: int = None, eps: float = None):
    for filename in tqdm(os.listdir(preds_dir)):
        try:
            uniprot_id = filename.split(".")[0]
            if not os.path.isfile(f"{distance_maps_dir}/{uniprot_id}.npy"):
                continue

            small_probs, metal_probs, nuclear_probs = parse_preds_file(os.path.join(preds_dir, filename))
        except:
            print(f'error parsing binding probabilities of file: {filename}')
            continue

        distance_map = np.load(f"{args.distancemap}/{uniprot_id}.npy")
        distances = distance_map[:, :, layer_index]
        valid_layer = distance_map[:, :, 4]

        small_refine = refine(small_probs, distances, valid_layer, epsilon=eps, k=k, l=1, overlap=True)
        metal_refine = refine(metal_probs, distances, valid_layer, epsilon=eps, k=k, l=1, overlap=True)
        nuclear_refine = refine(nuclear_probs, distances, valid_layer, epsilon=eps, k=k, l=1, overlap=True)

        write_to_file(filename=f"{uniprot_id + '.bindRefine_eps_' + str(eps) if k is None else uniprot_id + '.bindRefine_k_' + str(k)}", outdir=outdir,
                      small_hotspot=small_refine, metal_hotspot=metal_refine, nuclear_hotspot=nuclear_refine,
                      small_probs=small_probs, metal_probs=metal_probs, nuclear_probs=nuclear_probs)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This tool identifies a section of a protein with the highest average ligand binding probability.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required.add_argument('-p', '--predsdir', required=True,
                          help='directory containing predictions in specific format, see sample file. If your predictions are not available in this specific format. Please use the functions directly.')
    required.add_argument('-o', '--outdir', required=True, help='output directory')
    required.add_argument('-d', '--distancemap', required=True, help='directory containing protein distance maps, see sample file of distance map for required file structure and name.')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-eps', '--epsilon', type=float, help='set this value if you want to use epsilon mode of bindRefine')
    group.add_argument('-k', '--k', type=int, help='set this value if you want to use k mode of bindRefine')

    optional.add_argument('-l', '--layer', required=False, type=int, default=3, help='index of distance map layer. Options: 0 N, 1 C-alpha, 2 C-beta and 3 backbone C distances')

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

    main(outdir=args.outdir, preds_dir=args.predsdir, distance_maps_dir=args.distancemap, layer_index=args.layer, k=args.k, eps=args.epsilon)
