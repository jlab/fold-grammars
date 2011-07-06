== 1. Introduction ==
This repository contains all necessary files to build various RNA folding 
programs with the Algebraic Dynamic Programming (ADP) compiler "Bellmans GAP", 
written by Georg Sauthoff. Currently, there are four different Models / Grammars
for RNA secondary structures available:
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
    Currently, we provide two different signatures (subdir Signatures); their 
	prefix is "sig_". One is for grammars as used in RNAfold - namely NoDangle, 
	OverDangle and MicroState, thus this signature is called "sig_rnafold.gap". 
	The second is for the more complex grammar of RNAshapes - namely MacroState, 
	thus we call this signature "sig_rnashapes.gap".
 2. Algebra (different scoring or classification schemas for evaluating and 
    selecting from the search space)
    Prefix for algebras is "alg_". At this moment, we have implemented four 
	kinds of algebras:
    A. alg_rnafold_dotBracket: to build the Vienna Dot Bracket String for 
	   secondary structures
    B. alg_rnafold_mfe: to compute the structure of minimal free energy
    C. alg_rnafold_pfunc: to compute the partition function value for an 
	   ensemble of secondary structures
    D. alg_rnafold_shapes: this file actually contains five algebras, one for 
	   each shape abstraction level
    There exists extra versions of the algebras for the RNAshapes signature. 
	These are necessary, since both types of grammars use different algebra 
	functions and for MFE and Pfunc we also need special data types. In addition, 
	we have an experimental alg_rnashapes_centers algebra, where we get less 
	abstract than traditional shapes, because we also note their helix positions.
 3. Grammar (spans the whole searchspace of the given problem in a combinatorial 
    way)
    The subdirectory "Grammars" contains four grammars for RNA folding. Prefix is 
    "gra_". Note that NoDangle and OverDangle are identical. Together with 
	MicroDangle, these three grammars use the RNAfold signature, while 
	MacroStates uses signature RNAshapes.
 4. Instances (the actual appliance of the three former things, where different 
	combinations of algebras solve different tasks)
    Instances are defined in the "main" file of a Bellmans GAP program, namely 
	in "nodangle.gap", "overdangle.gap", "microstate.gap" and "macrostate.gap". 
	They also include the other modules and contain some technical instructions, 
	like including special functionality via C++ header files.

== 3. Prerequisites ==
Prerequisites for compiling an Instance are:
 - the Bellmans GAP compiler: http://www.gapc.eu/
   Georg kannst Du hier aufschreiben wie man da dran kommt?
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
