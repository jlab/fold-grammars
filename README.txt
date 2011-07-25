== 1. Introduction ==
This repository contains all necessary files to build various RNA folding 
programs with the Algebraic Dynamic Programming (ADP) compiler "Bellmans GAP", 
written by Georg Sauthoff. Currently (19.7.2011), there are four different 
Models / Grammars for RNA secondary structures available:
 - NoDangle
 - OverDangle
 - MicroState
 - MacroState
Their differences are discussed in detail in the 2011 BMC Bioinformatics Paper 
"Lost in folding space? Comparing four variants of the thermodynamic model for 
RNA secondary structure prediction." by Stefan Janssen et al.. Different tasks 
can be computed by different "algebras" or products of them, called 
"instances". 
For example:
 - prediction of a single, "optimal" MFE structure, i.e. its energy value and 
   its Vienna Dot Bracket representation: 
     mfe * dotBracket
 - computation of representative structures of different abstract shapes: 
     shape5 * mfe * dotBracket
 - computation of Boltzmann probabilities accumulated over all structures 
   of the same abstract shape: 
     shape5 * pfunc

== 2. Structure of repository ==
The source code is split into several modules. We use the ".gap" file extension. 
A Bellmans GAP program consists of four components:
 1. Signature (kind of an interface for grammars and algebras)
    Currently, we provide one signature (subdir Signatures); its prefix is 
	"sig_". "sig_foldrna.gap" fits for all four models, but not all functions, 
    defined in this signature, are really used in each model. For example, all 
	functions for the very complex MacroState grammar are defined in the 
	signature, but when we use NoDangle many of them remain unused. While this 
	happens on purpose it triggers some "Signature symbol xxx unused in grammar 
	yyy". Please, for this case, ignore them.
 2. Algebra (different scoring or classification schemas for evaluating and 
    selecting from the search space)
    Prefix for algebras is "alg_". At this moment, we have implemented four 
	kinds of algebras:
    A. alg_dotBracket: to build the Vienna Dot Bracket String for secondary 
	   structures
    B. alg_mfe: to compute the structure of minimal free energy
    C. alg_pfunc: to compute the partition function value for an ensemble of 
	   secondary structures
    D. alg_shapes: this file actually contains five algebras, one for each shape 
	   abstraction level
    There exists extra versions of the MFE and Pfunc algebras for the MacroState 
	model (alg_mfe_macrostate and alg_pfunc_macrostate). These extra versions 
	are necessary, since both algebras use special data types in the MacroState 
	model to keep track of different results for different dangling variants.
	
	E. Helix-Center (highly experimental)
	In addition, we have an experimental alg_centers algebra, where we get less 
	abstract than traditional shapes, because we also note their helix positions.
	A collection of alternative center definitions is provided by Jiabin Huang 
	(Uni Freiburg) in alg_hishapes.gap. To calculate the probability of a helix
	at a explicit position, we need a grammar that ensures that each search space
	candidate contains exactly one tagged helix (besides the completely unpaired
	structure). Similar to base pair probabilities, we want to know the probability
	of this helix independent of all other structural features of the remaining 
	bases for the candidates. That's the purpose of gra_macrostate_centers.gap
	in combination with macrostateCenter.gap.
 3. Grammar (spans the whole searchspace of the given problem in a combinatorial 
    way)
    The subdirectory "Grammars" contains four grammars for RNA folding. Prefix is 
    "gra_". Note that NoDangle and OverDangle are identical. 
 4. Instances (the actual appliance of the three former things, where different 
	combinations of algebras solve different tasks)
    Instances are defined in the "main" file of a Bellmans GAP program, namely 
	in "nodangle.gap", "overdangle.gap", "microstate.gap" and "macrostate.gap". 
	They also include the other modules and contain some technical instructions, 
	like including special functionality via C++ header files.
 5. The main directory contains some more misc files, namely 
    - pfunc_filter_foldrna.hh: augments the core Bellmans GAP system to 
	  heuristically filter low probability shape classes for the case of a 
	  (shape * pfunc) product.
	- pfunc_filter_macrostate.hh: same as pfunc_filter_foldrna.hh, but for the 
	  special data types for MacroStates.
	- pfunc_answer_macrostate.hh: special data types for MacroState algebras MFE 
	  and Pfunc.
    - sampling_wrapper.pl: a tiny perl wrapper for stochastic backtracing of 
	  shapes. Input is an RNA sequence, the number of samples and the shape level.
	  It uses preliminary compiled GAPc binaries to execute the sampling and 
	  returns the frequency of the shapes.
    - structure2shape.gap: a GAPc program, which parses a Vienna Dot Bracket 
	  string and returns a shape string of level 1 to 5. (not very efficient 
	  parsing)
    - README.txt: the file you are currently reading.
	
== 3. Prerequisites ==
Prerequisites for compiling an Instance are:
 - the Bellmans GAP compiler: http://www.gapc.eu/
 - Since Bellmans GAP compiler produces C++ target code, you need an C++ compiler. 
   We recommend "g++" and "make"
 - some editor to optionally modify the program.


== 4. How to compile an instance ==
Let us assume you want to build a program, which computes the best minimal free 
energy secondary structure, i.e. its energy value in kcal/mol plus its Vienna 
Dot Bracket representation for the NoDangle grammar:
 1. enter any console and cd to our repository
 2. generate C++ target code with Bellmans GAP compiler for the predefined 
    instance "mfepp" by typing: 
	  "gapc -i mfepp -t nodangle.gap"
 3. compile the generated target code via 
      "make -f out.mf"
 4. use the new binary 
      "./out ACUGACUAGCUAGUCGUUAC"
	  
See http://www.gapc.eu/faq.html for more information about Bellmans GAP compiler.
