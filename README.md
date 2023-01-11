# PyGMET
Python version of the Gridded Meteorological Ensemble Tool (GMET: https://github.com/NCAR/GMET). PyGMET has some technical differences compared to Fortran-based GMET.

# What does PyGMET do?
Note: All outputs are provided as netcdf4 files  
- Deterministic estimates of any input variable (defined in the configure file)
- Probabilistic estimates of any number of ensemble members
- Intermediate outputs: spatiotemporally correlated random fields, spatial correlation, temporal auto correlation (lag 1), nearby station index/weights, etc ...

# How to use PyGEMT?
- Set up Python environment  
Required packages should be introduced here in the form of requirements.txt or environment.yml

- Prepare the configuration file  
Use testcase.config.toml in ./src as the template to set up new cases

- Run PyGMET  
python main.py testcase.config.toml

# Existing PyGMET Datasets
EMDNA: Ensemble Meteorological Dataset for North America, https://doi.org/10.20383/101.0275  
EM-Earth: The Ensemble Meteorological Dataset for Planet Earth, https://doi.org/10.20383/102.0547

# Notes
This code is a work in progress and is provided without guarantee of fitness for any particular application.  
The PyGMET_SHARP branch is the most recent. Branch structure may be changed during the development.  

# References:
## PyGMET
Tang, G., Clark, M. P., & Papalexiou, S. M. (2022). EM-Earth: The Ensemble Meteorological Dataset for Planet Earth. Bulletin of the American Meteorological Society, 103(4), E996–E1018. https://doi.org/10.1175/BAMS-D-21-0106.1  
Tang, G., Clark, M. P., Papalexiou, S. M., Newman, A. J., Wood, A. W., Brunet, D., & Whitfield, P. H. (2021). EMDNA: an Ensemble Meteorological Dataset for North America. Earth System Science Data, 13(7), 3337–3362. https://doi.org/10.5194/essd-13-3337-2021
## Original GMET
Newman, A. J. et al. (2020) ‘Probabilistic Spatial Meteorological Estimates for Alaska and the Yukon’, Journal of Geophysical Research: Atmospheres, 125(22), pp. 1–21. doi: 10.1029/2020JD032696.   
Newman, A. J. et al. (2019) ‘Use of daily station observations to produce high-resolution gridded probabilistic precipitation and temperature time series for the Hawaiian Islands’, Journal of Hydrometeorology, 20(3), pp. 509–529. doi: 10.1175/JHM-D-18-0113.1.    
Newman, AJ, MP Clark, J Craig, B Nijssen, AW Wood, E Gutmann, N Mizukami, L Brekke, and JR Arnold, 2015, Gridded Ensemble Precipitation and Temperature Estimates for the Contiguous United States, J. Hydromet., doi: http://dx.doi.org/10.1175/JHM-D-15-0026.1  
Clark, M. P. and Slater, A. G. (2006) ‘Probabilistic Quantitative Precipitation Estimation in Complex Terrain’, Hydrometeorology, Journal O F, (2000), pp. 3–22.
