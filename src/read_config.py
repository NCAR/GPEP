# read model configurations

import toml, pathlib, json

def read_config(config_file):

    config = toml.load(config_file)

    # load model settings
    modelsettings_file = config['modelsettings_file']
    modelsettings_file = str(pathlib.Path(config_file).parent / modelsettings_file)
    settings = toml.load(modelsettings_file)

    config.update(settings)

    print('#'*50)
    print('Configuration file:', config_file)
    print(json.dumps(config, sort_keys=True, indent=4))
    print('#'*50, '\n'*2)

    return config