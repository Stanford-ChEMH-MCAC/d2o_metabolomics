# d2o_metabolomics

This repository contains R and python code for analyzing metabolomics experiments with heavy isotope (especially deuterium) labels.  The repo contains:
* `get_MID.py`, a command-line python script for extracting mass isotopologue distributions from `*.mzML` files, 
* several [jupyter notebooks](http://jupyter.org/install.html) with R code for analysis of MID data 
* a jupyter notebook for analysis of [`xcms`](https://bioconductor.org/packages/release/bioc/html/xcms.html) comparisons of labeled vs. non-labeled metabolomes analyzed by LC-MS.  

We are publishing the code not because it is a masterwork of software engineering, but because (a) we want to be transparent about how we analyze data that we report and publish on in the scientific literature, and (b) we hope that the code may be useful to someone else.

## Installation

No code here is packaged for easy installation, you will have to download the source and manually integrate it into your own stack.

## Requirements

Requirements are many and varied.  A rough guide:

* for `get_MID.py`:
	* python 3.5 or higher
	* [rdkit](http://www.rdkit.org/), 
	* [pyteomics](https://pythonhosted.org/pyteomics/) for reading `*.mzML` files
	* numpy, pandas, and other standard python libraries.
	
* for the R-based jupyter notebooks:
	* [xcms](https://bioconductor.org/packages/release/bioc/html/xcms.html)
	* [tidyverse](https://www.tidyverse.org/)
	* [viridis](https://cran.r-project.org/web/packages/viridis/index.html)
	* [janitor](https://cran.r-project.org/web/packages/janitor/index.html)
	* [ecipex](https://cran.r-project.org/web/packages/ecipex/index.html)

* general tip: use the [conda](https://conda.io/docs/) system for package management and installation.
	* installing `xcms` is probably an exception to above.

## Contributing

Pull requests are welcome.  If you use the code, please cite our paper (details forthcoming).  If you have questions about how any bit of this works, please file an [issue](https://github.com/Stanford-ChEMH-MCAC/d2o_metabolomics/issues) here on GitHub.

## Authors

* The code, as of May 2018, was written by [Curt R. Fischer](https://github.com/Stanford-ChEMH-MCAC/).

## License

This code is licensed with the [MIT License](https://opensource.org/licenses/MIT).
