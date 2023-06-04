# PyGMET

PyGMET is the Python version of the Gridded Meteorological Ensemble Tool ([GMET](https://github.com/NCAR/GMET)), a tool for generating gridded estimates for any meteorological variable using diverse methods. PyGMET comes with many methodological and technical differences and improvements compared to GMET.  

**Note: PyGMET is under active development. A formal release is expected to happen in late June.**  

<p align="center">
  <img src="https://github.com/guoqiang-tang/public_pictures/blob/main/california2017_ensemble_prcp.gif" alt="ensemble precipitation" />
</p>

## Functionality

PyGMET can perform the following tasks and output all results as netcdf4 files: 

-   Generate deterministic estimates of any input variable defined in the configuration file using various regression methods, including machine learning methods supported by scikit-learn.
-   Generate cross-validation outputs with evaluation results at station points.
-   Generate probabilistic estimates of any number of ensemble members.
-   Provide intermediate outputs, such as spatiotemporally correlated random fields, spatial correlation, temporal autocorrelation (lag 1), nearby station index/weights, and more.

## Installation

PyGMET is built on common Python packages such as SciPy, scikit-learn, and xarray. It can be run as long as the packages listed in environment.yml or requirements.txt are installed. You can create virtual environments using the following instructions: 
- Pip
```  
cd /your/path/of/PyGMET  
virtualenv PyGMET-env  
source PyGMET-env/bin/activate  
pip install -r requirements.txt  
```  
- Conda
```  
conda env create -f environment.yml  
conda activate PyGMET-env  
```

## Usage

To run PyGMET, follow these steps:

1.  Prepare the configuration file.  
Use configuration files in the `./test_cases` folder as templates. Refer to `How_to_create_config_files.md` for more details.
2.  Run PyGMET  
`python main.py /your/path/config_filename.toml`.
3. Batch run / operational run  
When producing a dataset in a target domain or testing different method choices, PyGMET can be run in many batches (e.g., month by month). We recommend users run PyGMET in two steps to improve efficiency:
    -   **Test run**: Run PyGMET on a test period or the first batch. Basic outputs (e.g., nearby station information, weights, and spatial correlation structures) will be generated and saved, which can be used by following batch runs.
    -   **Batch run**: Run PyGMET without changing `outpath_parent` so that PyGMET can find the outputs generated in the test run.

## Test Case

The test case for raw GMET is directly used here. To get the test case and run it, use the following commands:
```  
cd ../src  
python main.py ../test_cases/testcase.config.static.toml  
```
Jupyter Notebooks in the `./docs` folder can be used to visualize PyGMET test case outputs.

## Existing Public PyGMET Datasets

PyGMET has produced two publicly available datasets:

-   EMDNA: Ensemble Meteorological Dataset for North America, [https://doi.org/10.20383/101.0275](https://doi.org/10.20383/101.0275)
-   EM-Earth: The Ensemble Meteorological Dataset for Planet Earth, [https://doi.org/10.20383/102.0547](https://doi.org/10.20383/102.0547)  

## References:  
- PyGMET  
Tang, G., Clark, M. P., & Papalexiou, S. M. (2022). EM-Earth: The ensemble meteorological dataset for planet Earth. _Bulletin of the American Meteorological Society_, _103_(4), E996-E1018.  
Tang, G., Clark, M. P., Papalexiou, S. M., Newman, A. J., Wood, A. W., Brunet, D., & Whitfield, P. H. (2021). EMDNA: An ensemble meteorological dataset for North America. _Earth System Science Data_, _13_(7), 3337-3362.

- Original GMET  
Newman, A. J., Clark, M. P., Wood, A. W., & Arnold, J. R. (2020). Probabilistic spatial meteorological estimates for Alaska and the Yukon. _Journal of Geophysical Research: Atmospheres_, _125_(22), e2020JD032696.   
Newman, A. J., Clark, M. P., Longman, R. J., Gilleland, E., Giambelluca, T. W., & Arnold, J. R. (2019). Use of daily station observations to produce high-resolution gridded probabilistic precipitation and temperature time series for the Hawaiian Islands. _Journal of Hydrometeorology_, _20_(3), 509-529.  
Newman, A. J., Clark, M. P., Craig, J., Nijssen, B., Wood, A., Gutmann, E., ... & Arnold, J. R. (2015). Gridded ensemble precipitation and temperature estimates for the contiguous United States. _Journal of Hydrometeorology_, _16_(6), 2481-2500.     
Clark, M. P., & Slater, A. G. (2006). Probabilistic quantitative precipitation estimation in complex terrain. _Journal of Hydrometeorology_, _7_(1), 3-22.  
