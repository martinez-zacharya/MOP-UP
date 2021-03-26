Version : 0.9.1


**SETUP**

1. This has only been run on Ubuntu, however it might still be able to run on other Linux distros and MacOS.

2. Input ```$ Git clone https://github.com/martinez-zacharya/BipartitePipeline.git```

3. You need to install various things including

  	1. If using Ubuntu, build-essential ```$ sudo apt-get install build-essential```
  	2. SiLiX (http://lbbe.univ-lyon1.fr/-SiLiX-?lang=en)
	3. Miniconda (https://docs.conda.io/en/latest/miniconda.html)

4. Then you can enter ```$ conda env create -f environment.yml``` to create the necessary conda environment


5. You need a protein fasta file with no spaces in the sequence names, as well as a delimiter that serves as a cutoff between genome name and gene name, where to the left of the delimiter is the genome and to the right is the gene. This is incredibly important, don't overlook this! In the example below, the delimiter is "=". Also, make sure the delimiter is a valid character that can be used in a filename, unlike "&".

	Ex. >GCA_010101_@Peptidase2
	    MEKLKALLTSKKFWTLVAAIVAALT


**RUNNING**

To run, first make sure that diamond, silix and the python scripts that come with the package are in your curret working directory. The required parameters for the pipeline are

	1. Name of the run
	2. The full path to your protein fasta file
	3. The delimiter that separates the cutoff between Genome and Gene, such as '@'
	4. The full path to your desired output directory
	
There are several optional arguments you can add, which include

	1. --miniden Minimum percent identity to accept blast hits for building families. Default is 35%
	2. --minover Minimum percent overlap to accept blast hits for building families. Default is 80%
	3. --cpu Number of threads you want Diamond to use
	4. --db Whether or not you want to use our curated Microviridae database
	5. --singleton Add this flag to remove singletons
	6. --iter How many infomap iterations to run
	7. --connect Add this flag to keep only proteins that make connections in the output
	
An example command

```$ python3 BipartitePipeline.py TestRun /stor/home/Proteomes.fasta = /stor/work/TestRunOutput --cpu 12```


**OUTPUT**

Here is a brief description of the Micropipe output

	1. AllTitularProteins.fasta - One 'representative' protein from each protein cluster in a single fasta file
	2. nameofrunForCytoscape.csv - A network file for viewing with a program such as cytoscape
	3. nameofrunCytoscapeHelper.csv - When added to Cytoscape with the network file, allows you to differentiate between the protein cluster and subgroup nodes
	4. nameofrunMaster.csv - The master file that lists all of the genomes in the Micropipe run and their respective network ID and subgroup
	5. ProteinFamilies - A directory filled with fasta files for each protein cluster
	6. SubgroupMembers - A directory that contains text files for each subgroup with the subgroup members listed
