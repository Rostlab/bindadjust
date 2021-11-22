# Repository of Bachelor's Thesis "Refinement of Binding Residue Predictions Using Distance Maps"

### Introduction
This is the repository of my bachelor's thesis "Refinement of Binding Residue Predictions Using Distance Maps" at RostLab, completed on November 15, 2021.
This repo not only contains the entire codebase but also some example input and output files. 
If your input files are in the exact same format as the sample files, you can use the scripts main methods.
Otherwise, it is advised to incorporate the implemented functions in your own code. 

The application constitutes 3 independent modules:

- **bindViz.py**  uses PyMol to visualize binding residues on the protein 3D structures. It outputs a .gif file of a rotating protein annotated with its binding residue.
- **bindAdjust.py** modifies the small, metal and nuclear binding probabilities of all residues based on the distances between them. 
- **bindRefine.py** identifies one or several sections of the protein with the highest average binding probability. 

###How To bindViz.py

```
usage: bindViz.py [-h] -p PREDSDIR -t TRUES -o OUTDIR [-lm LIGANDMAP]
                  [-lt LIGANDTYPE] [-c CUTOFF] [-res RESOLUTION] [-fpr FPR]
                  [-s]

This tool visualizes binding residue on the 3D structure of proteins. Make
sure PyMol is installed on your system.

required arguments:
  -p PREDSDIR, --predsdir PREDSDIR
                        directory containing predictions in specific format,
                        see sample file. If your predictions are not available
                        in this specific format. Please use the functions
                        directly. (default: None)
  -t TRUES, --trues TRUES
                        file containing known binding residues, see sample
                        file. If your true values are not available in this
                        specific format. Please use the functions directly.
                        (default: None)
  -o OUTDIR, --outdir OUTDIR
                        output directory for visualizations (default: None)

optional arguments:
  -lm LIGANDMAP, --ligandmap LIGANDMAP
                        ligand map, maps UniProt sequences to PBD structures.
                        If no map is provided, stored map will be used.
                        (default: None)
  -lt LIGANDTYPE, --ligandtype LIGANDTYPE
                        ligand type to analyse, options are: small, metal and
                        nuclear (default: small)
  -c CUTOFF, --cutoff CUTOFF
                        cutoff used when converting float binding
                        probabilities to binary predictions (default: 0.5)
  -res RESOLUTION, --resolution RESOLUTION
                        resolution of the render, number of pixels for height
                        and width, always a square (default: 500)
  -fpr FPR, --fpr FPR   frames per rotation - number of frames generated for
                        one full rotation of the protein, more frames lead to
                        longer render times (default: 12)
  -s, --spectrum        visualize probabilities as continous color spectrum
                        instead of binary predictions (default: False)
  ```
    
###Example call using example files:

```
python3 bindViz.py -o <outdir> -lt small -p files/example_input/predictions -t files/example_input/binding_residues_2.5_small.txt -res 1000 -fpr 10
```

###Required Packages For bindViz.py
It is important that PyMol is running and that the following packages are installed in the enviroment:
- imageio
- tqdm
- pymol

###Example Visualizations
####Comparing predicted to annotated residues
![](files/example_output/bindViz/1s2k.gif)

####Binding probabilities as a color spectrum
![](files/example_output/bindViz/3q4o_spectrum.gif)


##How To bindAdjust.py

```
usage: bindRefine.py [-h] -p PREDSDIR -o OUTDIR -d DISTANCEMAP
                     (-eps EPSILON | -k K) [-l LAYER]

This tool identifies a section of a protein with the highest average ligand
binding probability.

required arguments:
  -p PREDSDIR, --predsdir PREDSDIR
                        directory containing predictions in specific format,
                        see sample file. If your predictions are not available
                        in this specific format. Please use the functions
                        directly. (default: None)
  -o OUTDIR, --outdir OUTDIR
                        output directory (default: None)
  -d DISTANCEMAP, --distancemap DISTANCEMAP
                        directory containing protein distance maps, see sample
                        file of distance map for required file structure and
                        name. (default: None)
  -eps EPSILON, --epsilon EPSILON
                        set this value if you want to use epsilon mode of
                        bindRefine (default: None)
  -k K, --k K           set this value if you want to use k mode of bindRefine
                        (default: None)

optional arguments:
  -l LAYER, --layer LAYER
                        index of distance map layer. Options: 0 N, 1 C-alpha,
                        2 C-beta and 3 backbone C distances (default: 3)
```
###Example call using example files:

```
python3 bindAdjust.py -o <outdir> -p files/example_input/predictions -d files/example_input/distance_maps -C 15
```

###Required Packages For bindAdjust.py
- tqdm
- numpy


##How To bindRefine.py

```
usage: bindRefine.py [-h] -p PREDSDIR -o OUTDIR -d DISTANCEMAP
                     (-eps EPSILON | -k K) [-l LAYER]

This tool identifies a section of a protein with the highest average ligand
binding probability.

required arguments:
  -p PREDSDIR, --predsdir PREDSDIR
                        directory containing predictions in specific format,
                        see sample file. If your predictions are not available
                        in this specific format. Please use the functions
                        directly. (default: None)
  -o OUTDIR, --outdir OUTDIR
                        output directory (default: None)
  -d DISTANCEMAP, --distancemap DISTANCEMAP
                        directory containing protein distance maps, see sample
                        file of distance map for required file structure and
                        name. (default: None)
  -eps EPSILON, --epsilon EPSILON
                        set this value if you want to use epsilon mode of
                        bindRefine (default: None)
  -k K, --k K           set this value if you want to use k mode of bindRefine
                        (default: None)

optional arguments:
  -l LAYER, --layer LAYER
                        index of distance map layer. Options: 0 N, 1 C-alpha,
                        2 C-beta and 3 backbone C distances (default: 3)
```
###Example call using example files:
```
python3 bindRefine.py -o <outdir> -p files/example_input/predictions -d files/example_input/distance_maps -k 10
```

###Required Packages For bindRefine.py
- tqdm
- numpy



