import configparser
import argparse

def read_config(filename, section):
    config = configparser.ConfigParser()
    config.read(filename)
    
    parameters = {}
    
    if section in config:
        for key, value in config[section].items():
            parameters[key.lower()] = value  # Convert all keys to lowercase
    
    return parameters

def main():
    parser = argparse.ArgumentParser(description="Run MDA simulations and save subsets to FASTA")
    parser.add_argument("-c", type=str, default="MDAsim.ini", help="Path to config.ini file")
    args = parser.parse_args()
    config_params = read_config(args.c, "Simulation")
    print(config_params['theta'])
    theta_value = config_params['theta']
    theta_type = type(theta_value).__name__
    print(f"Type of Theta: {theta_type}")
    # Print the configuration parameters
    print("Configuration Parameters:")
    for key, value in config_params.items():
        print(f"{key}: {value}")
    for key, value in config_params.items():
        # Determine the type of each value
        value_type = type(value).__name__
        print(f"{key}: {value_type}")
  
if __name__ == "__main__":
    main()
