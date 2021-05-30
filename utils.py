import os.path
from parse_utils import *
from sklearn.cluster import DBSCAN
from pymol import cmd
import numpy as np

def calc_stats(preds, trues):
    TP = list(trues & preds)
    FP = list(preds - trues)
    FN = list(trues - preds)

    # ID TP FP TN FN Precision Recall F1 MCC
    if len(preds) > 0:
        precision = len(TP) / len(preds)
    else:
        precision = 0.0
    if len(FN) + len(TP) > 0:
        recall = len(TP) / (len(FN) + len(TP))
    else:
        recall = 0.0

    F1 = len(TP) / (len(TP) + 0.5 * (len(FP) + len(FN)))
    return len(TP), len(FP), len(FN), precision, recall, F1

    # if len(preds) > 0:
    #    print(f"Precision: {suffix} {str(round(len(TP)/len(preds),3))}")
    # if len(FN) + len(TP) > 0:
    #    print(f"Recall: {suffix} {str(round(len(TP) / (len(FN) + len(TP)), 3))}")

    # cast values to int and substract offset and cast back to str


def vis_binding(preds, trues, offset, pdb_id, uniprot_id, outdir, suffix, fps):

    TP = list(trues & preds)
    FP = list(preds - trues)
    FN = list(trues - preds)

    if len(preds) > 0:
        print(f"Precision: {suffix} {str(round(len(TP) / len(preds), 3))}")
    if len(FN) + len(TP) > 0:
        print(f"Recall: {suffix} {str(round(len(TP) / (len(FN) + len(TP)), 3))}")

    # cast values to int and substract offset and cast back to str
    TP = [str(int(i) + offset) for i in TP]
    FP = [str(int(i) + offset) for i in FP]
    FN = [str(int(i) + offset) for i in FN]

    # print(f"TP {TP}")
    # print(f"FP {FP}")
    # print(f"FN {FN}")

    # get protein
    selection_name = "prot"

    cmd.fetch(pdb_id, name=selection_name, quiet=1)
    cmd.bg_color(color = "white")
    cmd.hide("everything")
    # cmd.show("prot")
    # show whole protein as grey surface
    cmd.show(selection=selection_name, representation="surface")
    cmd.color(color="grey", selection=selection_name)

    # color TP, FP and FN if they exist
    if len(TP) > 0:
        cmd.select(name="TP", selection="resi " + "+".join(TP))
        cmd.color(selection="TP", color="green")

    if len(FP) > 0:
        cmd.select(name="FP", selection="resi " + "+".join(FP))
        cmd.color(selection="FP", color="red")

    if len(FN) > 0:
        cmd.select(name="FN", selection="resi " + "+".join(FN))
        cmd.color(selection="FN", color="yellow")

    for i in range(int(360 / int(fps))):
        cmd.turn("y", fps)
        cmd.png(os.path.join(outdir, uniprot_id + "_" + pdb_id + "_" + str(i).zfill(3) + ".png"), ray=1, quiet=1)

    cmd.delete("*")

    images = []
    for file in sorted(os.listdir(outdir)):
        if not file.endswith(".png"):
            continue
        filename = os.path.join(os.curdir, outdir, file)
        # print(filename)
        images.append(imageio.imread(filename))
        os.remove(filename)

    imageio.mimsave(os.path.join(outdir, uniprot_id + "_" + pdb_id + suffix + ".gif"), images)

    for file in os.listdir():
        if file.endswith(".cif"): os.remove(file)

def filter_preds(distancemaps, uniprot_id, preds, eps=10, min_samples=3):

    distance_map = np.load(f"{distancemaps}/{uniprot_id}.npy")
    valid_layer = distance_map[:, :, 4]
    c_backbone_distances = distance_map[:, :, 3]


    # TODO: careful here, substracted one to convert from 1-based indices in preds to 0-based indices in nparray
    preds = sorted(list([int(i) - 1 for i in preds]))  # uniprot_to_preds['P23873']]))
    # sorted(list(preds))#uniprot_to_preds[sample_id]))
    # print(preds)

    # first, we throw out all preds we we dont have info in distance map

    # only use preds with valid distance
    valid_entries = [i for i in preds if valid_layer[i][i] != 0.0]

    # extract preds distances from matrix
    relevant_distances = c_backbone_distances[valid_entries, :][:, valid_entries]

    # new method to find min_samples
    # min_samples = int(np.sqrt(len(preds)))
    print(f"mean distance = {str(np.mean(relevant_distances))}")
    if np.mean(relevant_distances) < 7:
        print("thrown out")
        preds = [i + 1 for i in preds]
        return preds
    # if(eps == 0.0):
    #     preds = [i + 1 for i in preds]
    #     return preds


    # now we cluster the relevant distances
    # print(relevant_distances)
    # if the distances are invalid, we just return the preds
    if len(relevant_distances) == 0:
        # we have to add back plus 1
        preds = [i + 1 for i in preds]
        return preds
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(relevant_distances)
    labels = db.labels_

    # Number of clusters in labels, ignoring noise if present.

    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)

    # print(f"Noise has been filtered: {str(n_noise_/len(labels)).format(4)}")
    #dont do anything if too many clusters
  #  if n_clusters_ > 2:
  #      preds = [i + 1 for i in preds]
  #      return preds

    print('Estimated number of clusters: %d' % n_clusters_)
    print('Estimated number of noise points: %d' % n_noise_)

    # we remove the noise calculated by dbscan from our predictions
    # get indices of noise
    noise_indices = [i for i in range(len(labels)) if labels[i] == -1]

    # remove from valid entries
    newpred = np.delete(valid_entries, noise_indices)

    # convert back to 1-based for viz
    newpred = [i + 1 for i in newpred]
    return newpred