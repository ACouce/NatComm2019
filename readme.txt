COMPUTER SIMULATION:

The two files named "main_antimut.R" and "serial_culture_antimut.R" contain the basic code used to simulate how antimutator alleles invade a mutator, asexual population. All the relevant parameters, variables and processes are commented throughout the scripts. To run the model, execute the source file "main_antimut.R" in any R interpreter (e.g. R studio in Windows, Mac and Linux; bash console in Linux), such as:

> source('main_antimut.R')

The output is a figure with two panels. The first panel shows the invasion dynamics of a single antimutator allele though time for three different conditions. The second one shows the corresponding effective selection coefficients, as discussed in the text. Typical execution time is approximately 2-3 minutes on a 2.20 GHz Intel(R) Core(TM) i5-5200U CPU.


GENOME ANALYSES:

The file named "grantham.py" contains the basic code used to compute the Grantham score for all of the possible substitutions per codon in a bacterial genome. In this example, the code is applied to the genome of Escherichia coli K12 ("EcK12.fasta", retrieved from www.genoscope.cns.fr). The per codon Grantham scores are grouped to produce the average protein-disrupting effect of base substitutions elevated in 3 prominent types of mutators, as discussed in the text. To run the model, execute the source file "grantham.py" in any Python interpreter (e.g. bash console in Linux), such as:

> python grantham.py 

The output is the accumulated Grantham score and the number of sites targeted by the 3 considered spectra, per spectrum and per ORF. These data are displayed on the screen and also written to a file called "EcK12_results.dat". Typical execution time is approximately 30 seconds on a 2.20 GHz Intel(R) Core(TM) i5-5200U CPU.

REFERENCE:
Couce, A. & Tenaillon, O. (2019). Mutation bias and GC content shape antimutator invasions'. TBD.
