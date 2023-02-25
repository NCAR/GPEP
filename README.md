# PyGMET
PyGMET is the Python version of the Gridded Meteorological Ensemble Tool (GMET: https://github.com/NCAR/GMET). Although the Fortran version of GMET is used as the template, PyGMET has many technical differences and a flexible structure. Instead of focusing on precipitation and temperature, PyGMET can address any variable.  

# What does PyGMET docd 
- Deterministic estimates of any input variable (defined in the configure file)
- Probabilistic estimates of any number of ensemble members
- Intermediate outputs: spatiotemporally correlated random fields, spatial correlation, temporal auto correlation (lag 1), nearby station index/weights, etc ...  

All outputs are provided as netcdf4 files. 

# How to install PyGEMT
PyGMET does not have harsh requirement on Python environment. It can be run as long as the packages listed in environment.yml or requirements.txt are installed. You can also create virtual environments following below instructions. 
- pip    
<code>
cd /your/path/of/PyGMET  
virtualenv PyGMET-env  
source PyGMET-env/bin/activate  
pip install -r requirements.txt  
</code>

- conda    
<code>
conda env create -f environment.yml
conda activate PyGMET-env
</code>

# How to run PyGMET  
1. Prepare the configuration file  
Use testcase.config.toml in ./src as the template to set up new cases. The configuration file can be put anywhere.

2. Run the PyGMET  
python main.py config.toml  

# Test case for PyGMET  
The test case for raw GMET is directly used here. The below codes show how to get the test case and how to run the test case. 
<code>
cd ./test_cases  
python get_testcase.py  
cd ../src
python main.py ../test_cases/testcase.config.toml
</code>  

After running the test_cases, Jupyter Notebooks in the ./docs folder can be used to visualize PyGMET outputs.  

# Notes
This code is a work in progress and is provided without guarantee of fitness for any particular application.  
The PyGMET_SHARP branch is the most recent. Branch structure may be changed during the development.  

# Existing PyGMET Datasets
EMDNA: Ensemble Meteorological Dataset for North America, https://doi.org/10.20383/101.0275  
EM-Earth: The Ensemble Meteorological Dataset for Planet Earth, https://doi.org/10.20383/102.0547

# References:
## PyGMET
Tang, G., Clark, M. P., & Papalexiou, S. M. (2022). EM-Earth: The Ensemble Meteorological Dataset for Planet Earth. Bulletin of the American Meteorological Society, 103(4), E996–E1018. https://doi.org/10.1175/BAMS-D-21-0106.1  
Tang, G., Clark, M. P., Papalexiou, S. M., Newman, A. J., Wood, A. W., Brunet, D., & Whitfield, P. H. (2021). EMDNA: an Ensemble Meteorological Dataset for North America. Earth System Science Data, 13(7), 3337–3362. https://doi.org/10.5194/essd-13-3337-2021
## Original GMET
Newman, A. J. et al. (2020) ‘Probabilistic Spatial Meteorological Estimates for Alaska and the Yukon’, Journal of Geophysical Research: Atmospheres, 125(22), pp. 1–21. doi: 10.1029/2020JD032696.   
Newman, A. J. et al. (2019) ‘Use of daily station observations to produce high-resolution gridded probabilistic precipitation and temperature time series for the Hawaiian Islands’, Journal of Hydrometeorology, 20(3), pp. 509–529. doi: 10.1175/JHM-D-18-0113.1.    
Newman, AJ, MP Clark, J Craig, B Nijssen, AW Wood, E Gutmann, N Mizukami, L Brekke, and JR Arnold, 2015, Gridded Ensemble Precipitation and Temperature Estimates for the Contiguous United States, J. Hydromet., doi: http://dx.doi.org/10.1175/JHM-D-15-0026.1  
Clark, M. P. and Slater, A. G. (2006) ‘Probabilistic Quantitative Precipitation Estimation in Complex Terrain’, Hydrometeorology, Journal O F, (2000), pp. 3–22.
