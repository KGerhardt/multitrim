# multitrim

This is the development script for the new MiGA trimming approach.

To install the requirements, create a conda environment using multitrim.yml. Navigate to the directory in which multitrim.yml is place, and enter the following command:

conda env create -f multitrim.yml.

This will create a conda environment with the correct tools and will allow you to run the multitrim python script. Activate it using: 

conda activate multitrim

The python script can be run for paired end reads as: 

python3 multitrim.py -1 [FORWARD READS] -2 [REVERSE READS] --max -o [OUTPUT DIRECTORY]

For single-end reads, run it as: 

python3 multitrim.py -u [SE READS] --max -o [OUTPUT DIRECTORY]

NOTE:

Currently Falco is under development. The post-trim QC may fail. If this happens, the python script will issue an error and the post trim QC HTML will be missing. Please email me if this happens.

TEMPORARY:

I am in the process of creating an acutal workflow document for this script, but the following is a brief overview of the internal functioning of multitrim:

<hr />

#Workflow

* A subsample of up to 100K reads is taken from the input(s)
* The subsamples are run through FaQCs with report only mode on (no trimming) to detect adapters. Possible adapters come from this file: https://github.com/bio-miga/miga/blob/master/utils/adapters.fa
* The adapters detected (if any) are considered present if FaQCs reports them in a default 0.1% of reads. All adapters which are a part of the detected illumina kit(s) are included, e.g. detecting any one ilumina SE adapter will include ALL illumina SE adapters in the trim. The "families" of adapters can be seen at line breaks in the linked adapters.fa file
* Detected adapters are supplied to both FaQCs and Fastp in succession, so both tools attempt to trim adapters:
* First, FaQCs is run on the input reads with the -q 27 parameter, meaning that bases with < 27 quality score count against FaQCs' internal score parameter. This causes trimming to occur when enough <27 qual bases are found, and proceeds from both 5' and 3' ends separately.
* Second, Fastp is run on the post-FaQCs reads using a sliding window of size 3 and min avg. quality of 20. This is identical in behavior to trimmomatic's sliding window, but fastp is faster.
* Reads < 50 bp in length are removed.
* The final post-trim reads are output
* QC reports are performed on pre/post trim reads all at once.

<hr />

#Summary

* In essence: input reads -> FaQCs trims originals -> fastp  trims FaQCs outputs -> output reads

<hr />

#Trimmers

* Tools for trimming:

* FaQCs: https://github.com/LANL-Bioinformatics/FaQCs

* Fastp: https://github.com/OpenGene/fastp

#QC

* Falco: https://github.com/smithlabcode/falco/tree/master/src

#Sampling

* SeqTK: https://github.com/lh3/seqtk
