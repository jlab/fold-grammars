# fold-grammars
This repo is a fork of janssen labs fold-grammars, available at: https://github.com/jlab/fold-grammars
Collection of bgap code for RNA folding
For rudimentary documentation, consult https://bibiserv.cebitec.uni-bielefeld.de/fold-grammars/

# RNALoops

RNALoops introduced "motified" grammars and algebras, incorporating RNA 3D Motifs into RNA secondary structure predictions. The new algebras are available at /RNALoops/Algebras/Motifs/, the new grammars are in `/RNALoops/Grammars/`, being `gra_motified_macrostate.gap`, `gra_motified_microstate.gap`. The two grammars are slighly more complex versions of their respective base versions, introducing the option to identify RNA 3D Motifs surrounded only be a single base pair. Both grammars retain all qualities of their base versions `gra_macrostate` and `gra_microstate`.
Implementation of RNA 3D Motif recognition can be found under `/RNALoops/Extensions/motif.hh`, `/RNALoops/Extensions/MotifMap.hh`, and `/RNALoops/Extension/mot_header.hh`. The `mot_header.hh` file contains C style character arrays containing the currently available set of RNA 3D Motifs. 

Premade instances for our publication on `RNAmotiFold` are available in `/RNALoops/` under `RNAmotiFold.gap` and `RNAmotiFold_subopt_pfc.gap`. To replicate our results, use Bellman's GAP compiler arguments `--kbest --kbacktrace -t --no-coopt-class` when compiling the kbest MFE prediction algorithms. For subopt or partition function algorithms only use `-t`. All algorithms can also be installed through their associated MakeFiles, available under `/RNALoops/Misc/Applications/RNAmotiFold/RNAmotiFold` and `/RNALoops/Misc/Applications/RNAmotiFold/RNAmotiCes`. Simply call `make all -f makefile` while in the directory.

A python script for getting new versions of the RNA 3D Motif Atlas is available under `/RNALoops/Misc/Applications/RNAmotiFold/motifs/get_RNA3D_motifs.py`, which will also add the new motifs into the `mot_header.hh` file. After updating the motif sets it is necessary to (re)compile any RNAmotiFold algorithms, else they won't have the new motifs available. 

# Setup:</br>
+ If you are not familiar with Bellman's GAP and the Bellman's GAP compiler, I would recommend going to https://github.com/RNABioInfo/RNAmotiFold.git and have the python wrapper take care of any and all installation. Refer to the Readme there for more information.
+ Otherwise, a modified version of the original Bellman's GAP compiler (https://github.com/jlab/gapc) is available under https://github.com/RNABioInfo/gapcM.git. This can currently only be installed from source, using CMAke, or from my ubuntu launchpad PPA. Directions are given in the README there. It functions the same as the normal GAPC, except there are some alterations to the generated output and some changes that were necessary to implement RNA 3D Motifs (like extending the shape_alphabet to be able to incorporate Motifs into abstract RNA shapes). With the modified Bellman's GAP compiler it is also possible to incorporate RNA 3D Motifs into your own analysis through our classifiyng algebras!
</br>
If anything should not work for you when trying to implement RNALoops, RNAmotiFold or the gapcM, please feel free to reach out to me through my public e-mail.</br>
