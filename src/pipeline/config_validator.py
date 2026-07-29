from pathlib import Path
def check_config_options(config_data={}):
    # Check all options of config_data to ensure options are available/correct
    if len(config_data) == 0 or config_data is None:
        raise ValueError("[!] No configuration data provided.")
    
    # Check "mode"
    try:
        if config_data["mode"] == "paired-end-mode":
            pass
        elif config_data["mode"] == "merged-mode":
            pass
        elif config_data["mode"] == "unpaired-mode":
            pass
        else:
            raise ValueError("[!] Invalid mode specified. Use 'paired-end-mode', 'merged-mode', or 'unpaired-mode")
    except KeyError:
        raise ValueError("[!] No mode specified - Please specify a valid mode.")
    
    # Check "reference-fasta"
    try:
        with open(config_data["reference-fasta"],"r") as f:
            for line in f:
                if line.startswith(">"):
                    break
    except KeyError:
        raise ValueError("[!] No reference FASTA file specified - Please specify a valid reference FASTA file")
    
    except FileNotFoundError:
        raise ValueError("[!] Reference FASTA file not found")
    
    # Check "output-directory"
    try:
        config_data["output-directory"]
    except KeyError:
        raise ValueError("[!] No output directory specified - Please specify a valid output directory")

    # Check sequence files exist
    if config_data["mode"] == "paired-end-mode":
        try:
            with open(config_data["paired-end-mode"]["R1"], "r") as f:
                pass
        except FileNotFoundError:
            raise ValueError("[!] Read 1 file not found")
        
        try:
            with open(config_data["paired-end-mode"]["R2"], "r") as f:
                pass
        except FileNotFoundError:
            raise ValueError("[!] Read 2 file not found")
    elif config_data["mode"] == "merged-mode":
        try:
            with open(config_data["merged-mode"]["R1"], "r") as f:
                pass
        except FileNotFoundError:
            raise ValueError("[!] Merged file not found")
        
    # Check main analysis parameters are bool
    try:
        for key in config_data["core-parameters"]:
            if not isinstance(config_data["core-parameters"][key], bool):
                raise ValueError(f"[!] Invalid value for '{key}' - must be a boolean (True/False)")
    except KeyError:
        raise ValueError("[!] No core parameters specified")
    
    # Check processing parameters are valid, if processing parameters exist
    try:
        if config_data["core-parameters"]["do-processing"] == True:
            try:
                if not isinstance(config_data["processing-parameters"]["do-qtrimming"], bool):
                    raise ValueError("[!] Invalid value for 'do-qtrimming' - must be a boolean (True/False)")
                
                if config_data["processing-parameters"]["qtrimming-method"] not in ["rolling-trim", "simple-trim"]:
                    raise ValueError("[!] Invalid value for 'qtrimming-method' - must be 'rolling-trim' or 'simple-trim'")
                
                if not isinstance(config_data["processing-parameters"]["trimming-window-size"], int):
                    raise ValueError("[!] Invalid value for 'trimming-window-size' - must be an integer")
                
                if not isinstance(config_data["processing-parameters"]["trimming-quality-threshold"], int):
                    raise ValueError("[!] Invalid value for 'trimming-quality-threshold' - must be an integer")
            except KeyError:
                raise ValueError("[!] No processing parameters specified")
    except KeyError:
        # Non-bool "do-processing" caught by "core-parameters" check above, so this is only triggered if "do-processing" is not specified at all
        raise ValueError("[!] 'do-processing' parameter not specified in core parameters")