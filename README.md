Version : 0.6.0


**SETUP**

So this "pipeline" is still pretty ugly, hacky, and not terribly efficient, but it seems to work. There are only a few things you gotta do before you use it.

1. This has only been run on Ubuntu, however it might still be able to run on other Linux distros and MacOS.

2. You need to install various things including
  	1. If using Ubuntu (you should be), build-essential (sudo apt-get install build-essential)
  	2. SiLiX http://lbbe.univ-lyon1.fr/-SiLiX-?lang=en
  	3. Diamond http://www.diamondsearch.org/index.php?pages/installation/
  	4. Python 3 and several packages using pip/conda including
		1. Pip (sudo apt install python-pip)
		2. Pandas (pip3 install pandas)
		3. subprocess.run (pip3 install subprocess.run)
		4. InfoMap (pip3 install infomap)
		5. Figlet (pip3 install pyfiglet)
		6. BioPython (pip3 install biopython)

3. You need a protein fasta file with no spaces in the sequence names, as well as a delimiter that serves as a cutoff between genome name and gene name, where to the left of the delimiter is the genome and to the right is the gene. This is incredibly important, don't overlook this! In the example below, the delimiter is "=". Also, make sure the delimiter is a valid character that can be used in a filename, unlike "&".

	Ex. >GCA_010101_weight=1.5
	
	    >GCA_010101_weight=2.9


**RUNNING**

To run, first make sure that diamond is in the file you are currently in, as well as the python scripts that come with the package. The required parameters are
	1. Name of the run
	2. The full path to your protein fasta file
	3. The delimiter that separates the cutoff between Genome and Gene, such as =
	4. Type y or n depending on if you want to remove singleton protein clusters from the 		bipartite network
	5. Type y or n depending on if you want the output to only contain protein clusters that 	 connect to other protein clusters through a subgroup
	6. The full path to your desired output directory
	
There are several optional arguments you can add, which include
	1. --miniden Minimum percent identity to accept blast hits for building families. Default is	    35%
	2. --minover Minimum percent overlap to accept blast hits for building families. Default is 	    80%
	3. --cpu Number of threads you want Diamond to use
	
An example command

```$ python3 BipartitePipeline.py TestRun /stor/home/Proteomes.fasta = y n /stor/work/TestRunOutput --cpu 12```
