import argparse
import os

import imageio as imageio
import numpy as np
from tqdm import tqdm
from pymol import cmd

from utils import parse_trues, parse_preds, parse_ligand_map

def setup_pymol_load_structure(pdb : str):

    cmd.reinitialize()
    cmd.bg_color(color="white")

    # get protein structure
    selection_name = "prot"
    cmd.fetch(pdb, name=selection_name, quiet=1)

    # show only surface in grey
    cmd.hide("everything")
    cmd.show(selection=selection_name, representation="surface")
    cmd.set_color("TN", [239, 239, 239])
    cmd.color(color="TN", selection=selection_name)

    # define other colors
    cmd.set_color("TP", [49, 100, 255])
    cmd.set_color("FP", [254, 0, 0])
    cmd.set_color("FN", [136, 204, 238])

def render_and_save(outdir : str, pdb : str, resolution:int, fpr: int, filename_suffix : str = ''):
    images = []

    for i in range(fpr):
        cmd.turn("y", 360/fpr)
        cmd.png(os.path.join(outdir, pdb + "_" + str(i).zfill(3) + ".png"), ray=1, quiet=1, height=resolution, width=resolution)
        images.append(imageio.imread(os.path.join(outdir, pdb + "_" + str(i).zfill(3) + ".png")))
        os.remove(os.path.join(outdir, pdb + "_" + str(i).zfill(3) + ".png"))

    imageio.mimsave(os.path.join(outdir,  pdb + filename_suffix + ".gif"), images)

# if indices2 is passed, they are used as prediction and compared to indices which are considered trues
def vizResidues(pdb : str, outdir : str, indices, offset : int, indices2 = None, resolution : int = 500, fpr: int = 12):

    setup_pymol_load_structure(pdb)

    indices_with_offset = set([str(int(i) + offset + 1) for i in indices])

    if indices2 is None:
        cmd.select(name="True", selection="resi " + "+".join(indices_with_offset))
        cmd.color(selection="True", color="TP")
    else:
        indices2_with_offset = set([str(int(i) + offset + 1) for i in indices2])

        TP = list(indices_with_offset & indices2_with_offset)
        FP = list(indices2_with_offset - indices_with_offset)
        FN = list(indices_with_offset - indices2_with_offset)

        # color TP, FP and FN if they exist
        if len(TP) > 0:
            cmd.select(name="TP", selection="resi " + "+".join(TP))
            cmd.color(selection="TP", color="TP")

        if len(FP) > 0:
            cmd.select(name="FP", selection="resi " + "+".join(FP))
            cmd.color(selection="FP", color="FP")

        if len(FN) > 0:
            cmd.select(name="FN", selection="resi " + "+".join(FN))
            cmd.color(selection="FN", color="FN")

    # adjust camera and render mode
    cmd.set("ray_shadow", 0)
    cmd.set("ray_trace_mode", 1)
    cmd.move('z', -20)

    render_and_save(outdir, pdb, resolution, fpr)
    cmd.delete("*")
    if pdb + ".cif" in os.listdir(): os.remove(pdb + ".cif")

# list of available color spectra can be found here: https://pymolwiki.org/index.php/Spectrum
def vizSpectrum(pdb : str, outdir : str, probabilities: list(), offset : int, resolution : int = 500, fpr: int = 12, spectrum_color : str = "cyan_red"):

    setup_pymol_load_structure(pdb)

    # set all b values to 0
    cmd.alter("prot", "b=0")

    # set b values using probabilities
    for i, prob in enumerate(probabilities):
        color_value = probabilities[i]
        # this could be used to visualized folding: color_value = i / len(probabilities)
        cmd.alter("resi " + str(i + 1 + offset), "b = " + str(color_value))

    cmd.spectrum("b", spectrum_color, "prot")

    cmd.set("ray_opaque_background", 0)
    cmd.set("ray_shadow", 0)
    cmd.set("ray_trace_mode", 1)
    cmd.move('z', -20)

    render_and_save(outdir, pdb, resolution, fpr, filename_suffix="_spectrum")
    cmd.delete("*")
    if pdb + ".cif" in os.listdir(): os.remove(pdb + ".cif")


def main(outdir : str, preds_dir  : str, true_file : str, ligand_map : str, ligand_type : str, cutoff : float = 0.5, resolution: int = 500, fpr:int=10, mode: str = "binary"):

    print("parsing files...")
    uniprot_to_preds = parse_preds(filedir=preds_dir, ligand_type=ligand_type)
    uniprot_to_true = parse_trues(filepath=true_file)
    uniprot_to_pdb, uniprot_to_offset = parse_ligand_map(filepath=ligand_map)
    print("Done!")

    cmd.feedback("disable", "all", "everything")

    if mode == "binary":
        uniprot_ids = (uniprot_to_true.keys() & uniprot_to_preds.keys() & uniprot_to_pdb.keys())
        for uniprot_id in tqdm(uniprot_ids, desc=f"Rendering {len(uniprot_ids)} proteins in {mode = }", ncols=100):
            preds_indices = set([value[0] for value in np.argwhere(np.array(uniprot_to_preds[uniprot_id]) >= cutoff)])
            vizResidues(pdb = uniprot_to_pdb[uniprot_id], outdir=outdir, indices=uniprot_to_true[uniprot_id], offset=uniprot_to_offset[uniprot_id], indices2=preds_indices, resolution=resolution, fpr=fpr)

    if mode == "spectrum":
        uniprot_ids = (uniprot_to_preds.keys() & uniprot_to_pdb.keys())
        for uniprot_id in tqdm(uniprot_ids, desc=f"Rendering {len(uniprot_ids)} proteins in {mode = }", ncols=100):
            vizSpectrum(pdb=uniprot_to_pdb[uniprot_id], outdir=outdir, probabilities=uniprot_to_preds[uniprot_id], offset=uniprot_to_offset[uniprot_id], resolution=resolution, fpr=fpr)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This tool visualizes binding residue on the 3D structure of proteins. Make sure PyMol is installed on your system.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required.add_argument('-p', '--predsdir', required=True, help='directory containing predictions in specific format, see sample file. If your predictions are not available in this specific format. Please use the functions directly.')
    required.add_argument('-t', '--trues', required=True, help='file containing known binding residues, see sample file. If your true values are not available in this specific format. Please use the functions directly.')
    required.add_argument('-o', '--outdir', required=True, help='output directory for visualizations')

    optional.add_argument('-lm', '--ligandmap', required=False, help='ligand map, maps UniProt sequences to PBD structures. If no map is provided, stored map will be used.')
    optional.add_argument('-lt', '--ligandtype', required=False, default = "small", help='ligand type to analyse, options are: small, metal and nuclear')
    optional.add_argument('-c', '--cutoff', default=0.5, type=float, required=False, help='cutoff used when converting float binding probabilities to binary predictions')
    optional.add_argument('-res', '--resolution', default=500, type=int, required=False, help=('resolution of the render, number of pixels for height and width, always a square'))
    optional.add_argument('-fpr', '--fpr', default=12, type=int, required=False, help=('frames per rotation - number of frames generated for one full rotation of the protein, more frames lead to longer render times'))
    optional.add_argument('-s', '--spectrum', required=False, action='store_true', help="visualize probabilities as continous color spectrum instead of binary predictions")

    args = parser.parse_args()

    if args.spectrum:
        mode = "spectrum"
        print("visualizing probabilities as color spectrum. Cutoff is ignored.")
    else:
        mode = "binary"

    if not os.path.isdir(args.outdir):
        print('provided output directory does not exist or is not a directory.')
        exit(0)

    if not os.path.exists(args.trues):
        print('provided true file does not exist')
        exit(0)

    if not os.path.isdir(args.predsdir):
        print('provided prediciton directory does not exist or is not a directory.')
        exit(0)

    if args.ligandmap is None:
        print('no ligand map specified, using stored ligand map')
        ligand_map = "ligand_map_single_chain.txt"
    elif not os.path.exists(args.ligand_map):
        print('cannot find provided ligand map')
        exit(0)
    else:
        ligand_map = parser.ligand_map

    if args.ligandtype not in ['small', 'metal', 'nuclear']:
        print(f'{args.ligandtype} not a valid ligand type. ligand type has to be small, metal or nuclear.')
        exit(0)


    main(args.outdir, args.predsdir, args.trues, ligand_map, args.ligandtype, args.cutoff, resolution=args.resolution, fpr=args.fpr, mode=mode)