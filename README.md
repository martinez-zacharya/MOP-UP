Version : 0.9.0


**SETUP**

1. This has only been run on Ubuntu, however it might still be able to run on other Linux distros and MacOS.

2. You need to install various things including

  	1. If using Ubuntu, build-essential (sudo apt-get install build-essential)
  	2. SiLiX http://lbbe.univ-lyon1.fr/-SiLiX-?lang=en
  	3. Diamond http://www.diamondsearch.org/index.php?pages/installation/
  	4. Python 3 
	5. Several packages using pip (pip3 install -r requirements.txt)

3. You need a protein fasta file with no spaces in the sequence names, as well as a delimiter that serves as a cutoff between genome name and gene name, where to the left of the delimiter is the genome and to the right is the gene. This is incredibly important, don't overlook this! In the example below, the delimiter is "=". Also, make sure the delimiter is a valid character that can be used in a filename, unlike "&".

	Ex. >GCA_010101_=Peptidase2


**RUNNING**

To run, first make sure that diamond is in the file you are currently in, as well as the python scripts that come with the package. The required parameters are

	1. Name of the run
	2. The full path to your protein fasta file
	3. The delimiter that separates the cutoff between Genome and Gene, such as =
	4. The full path to your desired output directory
	
There are several optional arguments you can add, which include

	1. --miniden Minimum percent identity to accept blast hits for building families. Default is 35%
	2. --minover Minimum percent overlap to accept blast hits for building families. Default is 80%
	3. --cpu Number of threads you want Diamond to use
	4. --db Whether or not you want to use our curated Microviridae database
	5. --singleton Add this flag to remove singletons
	6. --iter How many infomap iterations to run
	7. --connect Add this flag to keep only proteins that make connections in the output
	8. --extra Add this flag to command MicroPipe to do the extra stuff
	
An example command

```$ python3 BipartitePipeline.py TestRun /stor/home/Proteomes.fasta = /stor/work/TestRunOutput --cpu 12```
