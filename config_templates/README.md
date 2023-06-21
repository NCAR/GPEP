# Configuration templates 

## testcase.config.dynamic.toml
This example uses ` ['latitude', 'longitude', 'elev', 'gradient_n_s', 'gradient_w_e']` as the predictors to perform locally weighted linear/logistic regression for the test case. To run the test case, following the codes as below:  
```
python main.py ../config_templates/testcase.config.dynamic.toml
```

## testcase.config.dynamic.toml
Same with testcase.config.dynamic.toml but also use precipitation and temperature as dynamic predictors in the regression

## testcase.config.RF.toml
This configuration file removes many optional parameters for simplicity. The regression method is Random Forest (RF). The configuration only performs spatial interpolation for precipitation, without activating ensemble estimation.

## model.settings.toml and model.settings.RF.toml
Model settings used in the configuration file. See ./docs/How_to_create_config_files.md for more detailed descriptions.