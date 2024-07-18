import configparser

def read_config(filename, section):
    config = configparser.ConfigParser()
    config.read(filename)
    
    parameters = {}
    
    if section in config:
        for key, value in config[section].items():
            parameters[key.lower()] = value  # Convert all keys to lowercase
    
    return parameters

def main():
    config_file = 'configs/MDAsim.ini'  # Path to your config file
    section = 'Simulation'

    config_params = read_config(config_file, section)
    
    # Print the configuration parameters
    print("Configuration Parameters:")
    for key, value in config_params.items():
        print(f"{key}: {value}")

if __name__ == "__main__":
    main()
