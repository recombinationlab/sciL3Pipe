
rule log_config:
    '''
    Copy config.yaml and place in logs folder with the date run
    '''
    output:
        out_dir + "logs/config_" + run_date + "yaml"
    run:
        with open(output[0], 'w') as out:
            yaml.dump(config, out, default_flow_style=False)