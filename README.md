# regdriver
regDriver to detect functional somatic mutaionts

## Required tools and modules

The following packages are required to run the regDriver method:
- Python 2.7: we recommend installing it from anaconda.com
- bedtools v25.9
- pybedtools (v0.7.8): for processing the annotation data files (e.g. conda install pybedtools=0.7.8 -c bioconda)
- The following python packages are also required to run the helper modules:
	- requests
	- userlib 
	- numpy
	

## Run the regDriver

	```python regDriver.py regDriver_params.conf```


	A copy of the regDriver_params.conf file that contains all parameters and options for the file can be found above 

