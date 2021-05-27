import os.path
import time
from utils import *
import pandas as pd
import argparse


def check_input(args):
    if not os.path.isdir(args.outdir):
        print("Output directory \"" + args.outdir + "\" does not exist.")
        exit(-1)
    if not os.path.isdir(args.predsdir):
        print("Output directory \"" + args.predsdir + "\" does not exist.")
        exit(-1)
    if not os.path.isfile(args.ligandmap):
        print("File \"" + args.ligandmap + "\" does not exist.")
        exit(-1)
    if not os.path.isfile(args.uniprot):
        print("File \"" + args.uniprot + "\" does not exist.")
        exit(-1)
    if not os.path.isfile(args.trues):
        print("File \"" + args.trues + "\" does not exist.")
        exit(-1)
    if args.distancemap is not None and not os.path.isdir(args.distancemap):
        print("Inout directory \"" + args.distancemap + "\" does not exist.")
        exit(-1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outdir', required=True, help='Output directory for files and visualizations')
    parser.add_argument('-u', '--uniprot', required=True, help='Input file: List of uniprot ids to analyse')
    parser.add_argument('-lm', '--ligandmap', required=True, help='Input file: Ligand map')
    parser.add_argument('-p', '--predsdir', required=True, help='Directory containing predictions')
    parser.add_argument('-t', '--trues', required=True, help='File containing true ligand indices and types')

    parser.add_argument('-th', '--threshold', required=False, help='Threshold used for classification of predictions', default=0.5)
    parser.add_argument('-fps', required=False, default=20, help='fps (not really fps) of rendered gif')

    parser.add_argument('-d', '--distancemap', required=False, help='Directory containing distance maps for clustering')
    parser.add_argument('-v', '--verbose', required=False, action='store_true', help="Verbose boolean")
    parser.add_argument('-s', '--save', required=False, action='store_true', help="Saves results as csv")
    parser.add_argument('-r', '--render', required=False, action='store_true', help="Renders gif of protein")
    # TODO: enum for small metal and nuclear? or pass column # of preds file

    # parse args
    args = parser.parse_args()

    # check if input files exist
    check_input(args)

    #booleans
    verbose = args.verbose
    render = args.render
    save = args.save

    #mandatory paths
    path_to_uniprots = args.uniprot
    path_to_ligand_map = args.ligandmap
    path_to_trues = args.trues
    path_to_preds_dir = args.predsdir
    outdir = args.outdir

    #optional path
    distance_maps_dir = args.distancemap

    #float and int
    threshold = args.threshold
    fps = args.fps

    if verbose:
        start_timer = time.time()
        print("Parsing input...", end="")

    # create empty data frame for statistics
    df = pd.DataFrame()

    # dictionaries for data
    # parse all input files
    uniprot_ids = parse_uniprot_ids(filepath=path_to_uniprots)
    uniprot_to_pdb, uniprot_to_offset = parse_ligand_map(filepath=path_to_ligand_map)
    uniprot_to_true = parse_trues(filepath=path_to_trues)
    uniprot_to_preds = parse_preds(filedir=path_to_preds_dir, threshold=threshold)

    if verbose:
        end_timer = time.time()
        print(" done after " + str(round(end_timer - start_timer, 2)) + " s.")

    ## write methode that prints stats :)
    # now vis the preds in pymol
    # we do this for every uniprot ID in the train set
    uniprot_ids = list(uniprot_ids)
    # epsilons = [10, 12, 15, 18, 20, 21, 22, 23, 24, 25, 26]
    # min_samples = [1, 2, 3, 4, 5, 6, 7]
    min_samples = [2]
    epsilons = [24]
    # cutoffs = [0.2,0.3,0.4,0.5]
    uniprot_ids = ['O43809', 'Q983T0', 'P82291']
    for min_sample in min_samples:
        for eps in epsilons:
            for uniprot_id in uniprot_ids:
                if uniprot_id not in uniprot_to_preds:
                    if verbose: print(f"Uniprot ID {uniprot_id} not in preds. Skipping...")
                    continue

                if not os.path.isfile(f"{distance_maps_dir}/{uniprot_id}.npy"):
                    if verbose: print(f"Uniprot ID {uniprot_id} has no distance map. Skipping...")
                    continue

                if uniprot_id not in uniprot_to_pdb:
                    if verbose: print(f"Uniprot ID {uniprot_id} not in ligand map. Skipping...")
                    continue

                if uniprot_id not in uniprot_to_true:
                    if verbose: print(f"Uniprot ID {uniprot_id} not in trues. Skipping...")
                    continue

                if render:
                    vis_binding(uniprot_id=uniprot_id, pdb_id=uniprot_to_pdb[uniprot_id], preds=uniprot_to_preds[uniprot_id], trues=uniprot_to_true[uniprot_id], offset=uniprot_to_offset[uniprot_id],
                                outdir=outdir, suffix="_before", fps=fps)

                TP, FP, FN, precision, recall, F1 = calc_stats(preds=uniprot_to_preds[uniprot_id], trues=uniprot_to_true[uniprot_id])

                if args.distancemap is not None:
                    preds_after_clust = filter_preds(distancemaps=args.distancemap, uniprot_id=uniprot_id, preds=uniprot_to_preds[uniprot_id], eps=eps, min_samples=min_sample)
                    preds_after_clust = set([str(i) for i in preds_after_clust])
                    TP2, FP2, FN2, precision2, recall2, F12 = calc_stats(preds=preds_after_clust, trues=uniprot_to_true[uniprot_id])
                    if render:
                        vis_binding(uniprot_id=uniprot_id, pdb_id=uniprot_to_pdb[uniprot_id], preds=preds_after_clust, trues=uniprot_to_true[uniprot_id],
                                    offset=uniprot_to_offset[uniprot_id],
                                    outdir=outdir, suffix="_after", fps=fps)
                    df = df.append([[uniprot_id, TP, FP, FN, precision, recall, F1, TP2, FP2, FN2, precision2, recall2, F12]])
                else:
                    df = df.append([[uniprot_id, TP, FP, FN, precision, recall, F1]])


            #stats = df.describe()
            #print(f"eps: {str(eps)} min_samples {min_sample} avg. precision: {stats[4][1]} avg. recall: {stats[5][1]} avg. F1: {stats[6][1]}")
            if save:
                print("saved .csv.")
                if args.distancemap is not None:
                    df.to_csv(f"{outdir}/small_eps{eps}_min_samp_{min_sample}.csv", sep="\t",
                              header=["ID", "TP_before", "FP_before", "FN_before", "Precision_before", "Recall_before", "F1_before", "TP", "FP", "FN", "Precision", "Recall", "F1"], index=False)
                else:
                    df.to_csv(f"{outdir}/small_without_cluster.csv", sep="\t", header=["ID", "TP", "FP", "FN", "Precision", "Recall", "F1"], index=False)
            df = pd.DataFrame()
