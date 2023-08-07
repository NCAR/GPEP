# read model configurations

import os, toml, pathlib, json

def read_config(config_file):

    # read main config file
    config = toml.load(config_file)

    # load model settings
    settings_file = config['modelsettings_file']
    settings_file = str(pathlib.Path(config_file).parent / settings_file)
    settings      = toml.load(settings_file)

    for k in ['flags', 'netcdf_info']:
        if k in settings:
            v = settings[k]
            settings.pop(k)
            settings.update(v)

    config.update(settings)

    # model paths
    for p in ['input_stn_list', 'input_stn_path', 'infile_grid_domain', 'outpath_parent', 'dynamic_predictor_filelist']:
        if p in config:
            pathconfig = str(pathlib.Path(os.path.abspath(config_file)).parent)
            if not os.path.isabs(config[p]):
                config[p] = pathconfig + '/' + config[p]

    print('#'*50)
    print('Configuration file:', config_file)
    
    if 'print_config' in config:
        if config['print_config'] == True:
            print(json.dumps(config, sort_keys=True, indent=4))
    print('#'*50, '\n'*2)

    return config
