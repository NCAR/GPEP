# read model configurations

import toml, pathlib, json

def read_config(config_file):

    # read main config file
    config = toml.load(config_file)

    # load model settings
    settings_file = config['settings_file']
    settings_file = str(pathlib.Path(config_file).parent / settings_file)
    settings      = toml.load(settings_file)

    config.update(settings)

    print('#'*50)
    print('Configuration file:', config_file)
    print(json.dumps(config, sort_keys=True, indent=4))
    print('#'*50, '\n'*2)

    return config