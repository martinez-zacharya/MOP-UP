Version : 0.1.1


**SETUP**

So this "pipeline" is still pretty ugly, hacky, and not terribly efficient, but it seems to work. There are only a few things you gotta do before you use it.

1. This has only been run on Ubuntu, however it might still be able to run on other Linux distros and MacOS.

2. You need to install various things including
  	1. If using Ubuntu (you should be), build-essential (sudo apt-get install build-essential)
  	2. Silix http://lbbe.univ-lyon1.fr/Download,3009.html?lang=fr
  	3. Diamond http://www.diamondsearch.org/index.php?pages/installation/
  	4. Python 3 and several packages using pip including
		1. Pip (sudo apt install python-pip)
		2. Pandas (pip3 install pandas)
		3. subprocess.run (pip3 install subprocess.run)
		4. InfoMap (pip3 install infomap)
		5. Figlet (pip3 install pyfiglet)
		6. PyInquierer (pip3 install PyInquirer)

3. You need a fasta file with no spaces in the sequence names, as well as a delimiter that serves as a cutoff between genome name and gene name, where to the left of the delimiter is the genome and to the right is the gene. This is incredibly important, don't overlook this! In the example below, the delimiter is "=".

	> Ex. GCA_010101_weight=1.5
	>     	GCA_010101_weight=2.9


**RUNNING**

To run, first make sure that diamond is in the file you are currently in, as well as the python scripts that come with the package. Then, just type

```$ python3 BipartitePipeline010.py```

and follow the instructions! Good luck :)
