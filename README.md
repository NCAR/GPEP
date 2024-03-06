# GPEP  

[![DOI](https://zenodo.org/badge/674783032.svg)](https://zenodo.org/badge/latestdoi/674783032)  

Geospatial Probabilistic Estimation Package (GPEP) is a python-based tool for generating gridded analyses of time-varying geophysical variables based on merging point/in-situ and spatially-distributed (i.e., gridded) observations and predictor variables. It was developed to expand on and advance the capabilities of the Gridded Meteorological Ensemble Tool (GMET: https://github.com/NCAR/GMET), which is written in FORTRAN. GPEP reproduces the baseline GMET capabilities (see Bunn et al, 2022) that were developed largely for meteorological dataset generation in climate and water resources applications, including creating inputs for hydrologic simulation and prediction. GPEP has a more flexible structure and provides a much broader array of methods than GMET, which relied solely on locally-weighted spatial linear and logistic regression methods. GPEP also has certain technical differences which were either unavoidable or pragmatic due to the conversion from FORTRAN to Python (including a different approach to cross-validation). GPEP development is currently supported by several US federal agency research programs. 

<p align="center">
  <img src="https://github.com/NCAR/GPEP/blob/develop/docs/california2017_ensemble_prcp.gif" alt="ensemble precipitation" />
</p>

## Functionality
GPEP can perform the following tasks: 
-   Generate deterministic estimates of any input variable defined in the configuration file using various regression methods, including machine learning methods supported by scikit-learn.
-   Generate cross-validation outputs with evaluation results at station points.
-   Generate probabilistic estimates of any number of ensemble members.
-   Provide intermediate outputs, such as spatiotemporally correlated random fields, spatial correlation, temporal autocorrelation (lag 1), nearby station index/weights, and more.

## Installation
GPEP attempts to rely on common library dependencies in Python (version 3 or greater) environment. It can be run as long as the packages listed in environment.yml or requirements.txt are installed. Users can also create virtual environments following below instructions. 
- Pip
```  
cd /your/path/of/GPEP  
virtualenv GPEP-env  
source GPEP-env/bin/activate  
pip install -r requirements.txt  
```  
- Conda
```  
conda env create -f environment.yml  
conda activate GPEP-env  
```

## Usage
Please refer to the demonstration notebook (`./docs/GPEP_demo.ipynb`), which includes downloading test cases, running the test cases, and visualizing the results. This notebook demonstrates a simple case of using GPEP.  

It is recommended to directly run GPEP using Python, follow these steps:

1.  Prepare the configuration file.  
Use configuration files in the `./config_templates` folder as templates. Refer to `./docs/How_to_create_config_files.md` for more details. Configuration files are the core of GPEP cases.  
2.  Run GPEP  
`python main.py /your/path/config_filename.toml`.
3. Batch run / operational run  
When producing a dataset in a target domain or testing different method choices, GPEP can be run in many batches (e.g., month by month). We recommend users run GPEP in two steps to improve efficiency:
    -   **Test run**: Run GPEP on a test period or the first batch. Basic outputs (e.g., nearby station information, weights, and spatial correlation structures) will be generated and saved, which can be used by following batch runs.
    -   **Batch run**: Run GPEP without changing `outpath_parent` so that GPEP can find the outputs generated in the test run.

## Test Cases for GPEP

The test cases of GPEP are open access on Zenodo: https://zenodo.org/record/8222852. They can be obtained using the `./tools/get_testcase.py` script, which is used in `./docs/GPEP_demo.ipynb`. See ./tools/README.md for more descriptions.  

## Notes
This code is a work in progress and is provided without guarantee of fitness for any particular application.  
The main branch is the formal released version and the develop branch is the most recent. Branch structure may be changed during the development.  

## GPEP Development Notes
GPEP development was initiated at NCAR to supercede GMET, a software created and developed at NCAR over multiple years within multiple projects supported primarily by Reclamation and USACE, two US federal water management agencies that were interested in improving the surface forcing datasets available for hydrological modeling, and particularly their characterization of uncertainty, which propagates into hydrologic forecast uncertainty.  

A 2021 USACE project (NCAR PI A. Wood) focused on characterizing the robustness of hydrological models for climate change analysis provided an opportunity to transition GMET from FORTRAN to Python in order to broaden the user base and ease of experimentation with new methods. Additional support was provided through the Reclamation Science and Technology (S&T) program within two projects oriented to develop surface forcings in the Big Thompson and Sacramento River basins. Development on the GMET code was halted.  

The python version of GMET, re-named to a more generalized 'GPEP', was subsequently created at NCAR and designed for backward compatibility with the usage modes and input configurations in the existing NCAR USACE and Reclamation GMET applications. GPEP includes both GMET regression and ensemble generation functionality while also expanding regression options as noted above. GPEP development relied on GMET test cases developed for both USACE and Reclamation projects at NCAR (led by A. Wood). The GPEP basic (i.e., GMET) regression functionality also incorporated some code from various python regression scripts previously written by G. Tang at the University of Saskatchewan, sponsored by the Global Water Futures Project. Those had been used in combination with the GMET (FORTRAN) ensemble generation code to generate several large-domain meteorological datasets (EMDNA and EM-Earth; see references below). 

In 2024-2025, GPEP development continues with support of projects from USACE, Reclamation and recently NOAA, which is sponsoring an effort in collaboration with U. of Oklahoma/NSSL to develop a GPEP-based 2-km US-wide ensemble multi-decadal to real-time surface meteorological forcing dataset.  

## Sponsorship
- (current) US Army Corps of Engineers, Climate Preparedness and Resilience Program (NCAR projects)
- (current) US Bureau of Reclamation, Snow Water Supply Forecasting Program; and Science and Technology Program (NCAR projects)  
- (current) NOAA OAR, Community Observations and Modeling (COM) program (NCAR projects)
- (prior, related) Global Water Futures (U. Saskatchewan project)

## Related References
### GPEP reference
- Tang, G, AW Wood, AJ Newman, MP Clark, and SM Papalexiou, 2023. GPEP v1.0: a Geospatial Probabilistic Estimation Package to support Earth Science applications. Geosci. Mod. Dev.  https://doi.org/10.5194/gmd-2023-172, 2024. 

### GMET (FORTRAN-based) datasets and methods
- __GMET v2.0__:  Bunn, PTW, AW Wood, AJ Newman, H Chang, CL Castro, MP Clark and JR Arnold, 2022, Improving station-based ensemble surface meteorological analyses using numerical weather prediction:  A case study of the Oroville Dam crisis precipitation event. J. Hydromet. 23(7), 1155-1169. https://doi.org/10.1175/JHM-D-21-0193.1  
- Liu, Hongli, AW Wood, AJ Newman and MP Clark, 2021, Ensemble dressing of meteorological fields: using spatial regression to estimate uncertainty in deterministic gridded meteorological datasets, AMS J. Hydromet., https://doi.org/10.1175/JHM-D-21-0176.1   
- Newman, A. J. et al. (2020) Probabilistic Spatial Meteorological Estimates for Alaska and the Yukon, Journal of Geophysical Research: Atmospheres, 125(22), pp. 1–21. https://doi.org/10.1029/2020JD032696  
- Newman, A. J. et al. (2019) Use of daily station observations to produce high-resolution gridded probabilistic precipitation and temperature time series for the Hawaiian Islands, Journal of Hydrometeorology, 20(3), pp. 509–529. https://doi.org/10.1175/JHM-D-18-0113.1  
- Mendoza, PA, AW Wood, EA Clark, E Rothwell, MP Clark, B Nijssen, LD Brekke, and JR Arnold, 2017, An intercomparison of approaches for improving predictability in operational seasonal streamflow forecasting, Hydrol. Earth Syst. Sci., 21, 3915–3935, 2017 (used a real-time implementation of GMET)
- Newman, AJ, MP Clark, J Craig, B Nijssen, AW Wood, E Gutmann, N Mizukami, L Brekke, and JR Arnold, 2015, Gridded Ensemble Precipitation and Temperature Estimates for the Contiguous United States, J. Hydromet., doi: http://dx.doi.org/10.1175/JHM-D-15-0026.1   
- Clark, M. P. and Slater, A. G. (2006) Probabilistic Quantitative Precipitation Estimation in Complex Terrain, Hydrometeorology, Journal O F, (2000), pp. 3–22. https://doi.org/10.1175/JHM474.1  

### Other datasets created using python regression scripts and GMET ensemble generation
- Tang, G., Clark, M. P., & Papalexiou, S. M. (2022). EM-Earth: The Ensemble Meteorological Dataset for Planet Earth. Bulletin of the American Meteorological Society, 103(4), E996–E1018. https://doi.org/10.1175/BAMS-D-21-0106.1  
- Tang, G., Clark, M. P., Papalexiou, S. M., Newman, A. J., Wood, A. W., Brunet, D., & Whitfield, P. H. (2021). EMDNA: an Ensemble Meteorological Dataset for North America. Earth System Science Data, 13(7), 3337–3362. https://doi.org/10.5194/essd-13-3337-2021

## Contacts
- Guoqiang Tang (guoqiang@ucar.edu), primary developer
- Andy Wood (andywood@ucar.edu), project(s) lead

## Additional Acknowledgements
We would like to acknowledge high-performance computing support from Cheyenne (doi:10.5065/D6RX99HX) and Derecho provided by NCAR's Computational and Information Systems Laboratory, sponsored by the National Science Foundation.






