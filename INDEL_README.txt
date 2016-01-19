==========================================================================================
Phasterate
==========================================================================================
PHAST: PHYLOGENETIC ANALYSIS WITH SPACE/TIME MODELS

EXTENDED WITH INDELS
==========================================================================================

This is an extended version of Phast, for any question please email:
"mira.han" atSymbol "unlv.edu"

Installation:
For installation instructions see: README.txt in this directory.

This serves as a quick start guide and reference. Please see PhastDocumentation.pdf
for more information.


Sample Runs:
Quick overview of sample data and running it. For further information please see
IndelDocumentation.pdf in this folder.

The following commands are executed from this directory PhastIndel-1.3/, we assume this
is the current working directory.

There are several sample file in SampleRuns with expected input and output files:

1) Running PhyloFit on a single file with the F84E model:

   bin/phyloFit --tree SampleRuns/SingleTest/ENSGT00390000000013.newick \
   SampleRuns/SingleTest/ENSGT00390000000013.fa --subst-mod F84E -G -O branches

   This will produce two output files:
   phyloFit.mod      //mod file containing the fitted model and tree.
   phyloFit.infoX     //infoX file containing the fitted rate matrix parameters.

   Copies of these files can be found in: SampleRuns/ModFilesData/SingleTest/

2) Running PhyloFit on multiple files with the F84E model:

   bin/phyloFit -M SampleRuns/DoubleNewickTest/ SampleRuns/DoubleMsaTest/ \
   --subst-mod F84E -G -O branches

   This will produce a folder phyloFitResults/ containing the following files:

   ENSGT00390000000013.mod  ENSGT00390000000013.mod.infoX
   ENSGT00390000000046.mod  ENSGT00390000000046.mod.infoX

   Copies of these files can be found in: SampleRuns/ModFilesData/DoubleTest/

Note: It is not recommended to use phyloP when the data was fitted with very little
      data. Our model was fitted with ~1000 files.

3) Running PhyloP to compute indel conservation scores for an alignemnt.

   bin/phyloP  --refidx 0 --wig-scores --mode CONACC --method LRT \
   SampleRuns/ModFilesData/DoubleTest/ENSGT00390000000013.mod \
   SampleRuns/SingleTest/ENSGT00390000000013.fa \
   -x SampleRuns/ModFilesData/DoubleTest/ENSGT00390000000013.mod.infoX

   Notice this is the file we fitted on step 2). The P-Scores are printed to standard
   output. A file computedLikelihoods.txt is also created containing the null, alt, and
   scales for each site in the alignment.

   Copies of these files can be found in: SampleRuns/OutputPhyloP/DoubleTest

4) Running PhyloP to compute indel conservation scores for an alignemnt using
   our ~1000 file fitted model:

   bin/phyloP  --refidx 0 --wig-scores --mode CONACC --method LRT \
   SampleRuns/ModFilesData/TotalData/ENSGT00390000000013.mod \
   SampleRuns/SingleTest/ENSGT00390000000013.fa \
   -x SampleRuns/ModFilesData/TotalData/infoX.txt

   Copies of these files can be found in: SampleRuns/OutputPhyloP/TotalData/

   Our Fitted rate matrix can be found in SampleRuns/ModifiedData/TotalData/
   To use with other alignments the user may simply copy and paste the rate matrix
   to a new properly names mod file with the tree for the alignment they wish to run.
   PhyloP actually gets the rate matrix information from the infoX file, not the mod
   file.
==========================================================================================
