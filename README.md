Version : 0.9.5


**SETUP**

1. This has only been run on Ubuntu, but it might work on other Linux distros

2. Input ```$ git clone https://github.com/martinez-zacharya/MicroPipe.git``` into the command line

3. You need to install various things including

  	1. If using Ubuntu, build-essential ```$ sudo apt-get install build-essential```
  	2. SiLiX (https://lbbe-web.univ-lyon1.fr/fr/SiLix)
		1. Make sure to place the silix file into the same directory as mop-up.py and make it exectuable
	3. Miniconda (https://docs.conda.io/en/latest/miniconda.html)
4. Enter the mop-up directory, enter ```$conda activate
5. Then you can enter ```$ conda env create -f environment.yml``` to create the necessary conda environment


6. You need a protein fasta file with 
	1. No spaces in the headers
	2. A delimiter that serves as a cutoff between genome name and gene name, where to the left of the delimiter is the genome and to the right is the gene. In the example below, the delimiter is "@". Also, make sure the delimiter is a valid character that can be used in a filename, unlike "&".

	> \>ClusterI_0_SCPE01000001_Alistipes_sp@gene_000001
	> MEKLKALLTSKKFWTLVAAIVAALTAFFTTSCTGYLKFRREGVHHDTVRYEQVIKHKNYSAWLSNQIDRSSWRRPMMLSVCSSGTAFLSRSSGVSISSYLISPPLIHSVLPSILSSSINSVSFITLNRNCLEQITPLFTSNLPKRKSLGVCLALMMIGVPLSRVSPSSLDERRRGGRKGRPRPKIRNIPIGGSHL
	
	3. No duplicated gene names



**RUNNING**

To run, first make sure that diamond, silix and the python scripts that come with the package are in your curret working directory. The required parameters for the pipeline are

	1. Name of the run
	2. The full path to your protein fasta file
	3. The delimiter that separates the cutoff between Genome and Gene, such as '@'
	4. The full path to your desired output directory
	
There are several optional arguments you can add, which include

	1. --miniden Minimum percent identity to accept blast hits for building families. Default is 35%
	2. --minover Minimum percent overlap to accept blast hits for building families. Default is 80%
	3. --cpu Number of threads you want Diamond to use. Default is all available threads
	4. --noMicro By adding this argument, you can elect to leave out the Microviridae database from the analysis
	5. --singleton Add this flag to remove singletons
	6. --iter How many infomap iterations to run. Default number is 1000
	7. --block-size Change block size for Diamond to increase performance depending on available RAM
	8. --sensitivity Tune alignment sensitivity for Diamond. Options include mid-sensitive, sensitive, more-sensitive, very-sensitive, and ultra-sensitive

An example command

```$ python3 mop-up.py TestRun /stor/home/Proteomes.fasta @ /stor/work/TestRunOutput --cpu 12```

**OUTPUT**

Here is a brief description of the Micropipe output

	1. nameofrunAllTitularProteins.fasta - One 'representative' protein from each protein cluster in a single fasta file
	2. nameofrunConnectingTitularProteins.fasta - One 'representative' protein from each protein cluster that connects to another cluster through a subgroup in a single fasta file
	3. nameofrunForCytoscape.csv - A network file for viewing with a program such as cytoscape
	4. nameofrunCytoscapeHelper.csv - When added to Cytoscape with the network file, allows you to differentiate between the protein cluster and subgroup nodes
	5. nameofrunMaster.csv - The master file that lists all of the genomes in the Micropipe run and their respective network ID and subgroup
	6. nameofrunProteinFamilies - A directory filled with fasta files for each protein cluster
	7. nameofrunSubgroupMembers - A directory that contains text files for each subgroup with the subgroup members listed
