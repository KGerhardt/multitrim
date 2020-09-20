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
