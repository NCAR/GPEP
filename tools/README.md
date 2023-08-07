# Tools for various purposes  
## Get test cases
The test cases of GPEP are open access on Zenodo: https://zenodo.org/record/8222852. Users can download them from the website, or simply use the script `get_testcases.py` provided here.  
`get_testcases.py` accepts two optional parameters: `-o <outpath>` where `<outpath>` can be replaced by the directory where test cases will be saved; `-n <x>` where `<x>` defines what test cases will be downloaded. For example, replacing `<x>` by `0,1` will download test cases 0 and 1, which correspond to the cali2017 (California prcp, tmean and trange) and UpCO_SWE (Upper Colorado snow water equivalent).    

Example usage - 1:  
`./tools/get_testcases.py`, which will download and uncompress the test case 0 to the default folder (i.e., ../GPEP_test_cases)

Example usage - 2:  
`./tools/get_testcases.py -o /your/path/GPEP_testcases -n 0,1`, which will download and uncompress the test cases 0 and 1 to the `/your/path/GPEP_testcases` folder  


## Other tools
- plot_ensemble_prcp.ipynb. This notebook shows how to plot a gif figure of ensemble precipitation
- plot_cross_validation_metrics.ipynb. This notebook shows how to plot cross validation results
