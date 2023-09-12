#!/bin/sh -l

#### The commands and scripts used were written by Carmen L. Wickware, Julian Trachsel, and Timothy A. Johnson

### Step 1: hmmsearch command in HMMER software was used to annotate the metagenome ###

# This command will require a library that is now available through NCBI
# The model and library (version 15) that were used were originally hosted on the JCVI site
# This information is now archived in NCBI should you want to use the exact version

echo "start hmmer" >> log
date >> log

hmmsearch -o BMD.TIGR --tblout BMD.TIGR.tsv --cut_tc TIGRFAMs_15.0_HMM.LIB metaspades_500_scaff.prot.fasta.ns

echo "finished hmmer" >> log
date >> logs

### Step 2: An annotation file (.gff) was created using python script TIGR2gff.py available in this repository (Py_R_Scripts) ###

# This will format the .tsv file created above into the annotation .gff file needed in subsequent steps
# The script can be run with the following commands

python TIGR2gff.py BMD.TIGR.tsv TIGR.gff

### Step 3: using featureCounts, part of the Subread software, features such as genes were counted ###

# This is done via a loop to create the command for each feature and run in parallel to speed up the process

echo "start featureCount-TIGR" >> log
date >> log

echo "start featureCount-TIGR" >> log
for i in ./mapping_files/*.bam; do echo "featureCounts -t func_gene -g gene_id -a TIGR.gff -o $i.counts $i"; done > featureCounts-command.sh
cat featureCounts-command.sh | parallel

echo "finished featureCount-TIGR" >> log
date >> logs

### Steps 4-5 remaining steps were done using R, with scripts available in this repository (Py_R_Scripts) ###

# Step 4 used the script TIGR_countmerge.R to combine all feature counts into one file

# Step 5 used the script TIGRfam_annotationMerge.R to combine the roles and subroles into each feature annotation
# This last step requires files that contain information on the broader roles of each TIGRfam annotation
# The files used in the Rscript are available in this repository (TIGRfam_Annotations_Roles)
