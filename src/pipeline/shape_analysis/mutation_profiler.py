import math
from pathlib import Path
import matplotlib.pyplot as plt

def data_retriever(mpileup_analysis:Path):
    # Data out dict
    data_out = {
        'positions':[],
        'ident_rates':[],
        'mut_rates':[],
        'ins_rates':[],
        'del_rates':[],
        'depth':[]
    }
    # Map keys to column positions in mpileup_analysis
    file_map = {
        "positions": 0,
        "ident_rates": 1,
        "depth": -1,
        "mut_rates": -2,
        "ins_rates": -4,
        "del_rates": -3,
    }
    # Open analysis file
    with open(mpileup_analysis, "r") as f:
        # Skip header
        next(f)
        for line in f:
            # Parse tsv fields
            fields = line.strip().split("\t")
            for key, index in file_map.items():
                data_out[key].append(float(fields[index]))

    return data_out

def overlay_plot_generator(data1:dict, data2:dict, output_directory):
    # Stack plots
    # Generate and save the rate_identical plot
    plt.plot(data1["positions"], data1["ident_rates"], color='red', label='modified')
    plt.plot(data2["positions"], data2["ident_rates"], color='blue', label='unmodified')
    plt.title('Rate Identical to Reference\nmodified v unmodified')
    plt.xlabel('Position')
    plt.ylabel('Rate Identical (identical base per read)')
    plt.xlim(0, max(data1["positions"]))
    plt.ylim(0, 1.2)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_directory}/Overlaid_Rate_identical_plot.png")
    plt.clf()
    
    # Generate and save the rate_mutation plot
    plt.plot(data1["positions"], data1["mut_rates"], color='red', label='modified')
    plt.plot(data2["positions"], data2["mut_rates"], color='blue', label='unmodified')
    plt.title('Rate Mutation from Reference')
    plt.xlabel('Position')
    plt.ylabel('Mutation rate (mutations per read)')
    plt.xlim(0, max(data1["positions"]))
    plt.ylim(0, 2.2) # The highest possible mutation rate / position is 2 (1 mismatch and 1 indel in every read at position)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_directory}/Overlaid_Rate_mutation_plot.png")
    plt.clf()

    # Generate and save the sequencing depth plot
    plt.plot(data1["positions"], data1["depth"], color='red', label='modified')
    plt.plot(data2["positions"], data2["depth"], color='blue', label='unmodified')
    plt.title('Sequencing Coverage Across Positions')
    plt.xlabel('Position')
    plt.ylabel('Depth')
    plt.xlim(0, max(data1["positions"]))
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_directory}/Overlaid_Depth_plot.png")
    plt.clf()

    # TODO - skipping indel overlay since its not very useful methinks

def mutation_profiler(m_full_analysis:Path, u_full_analysis:Path, output_directory):
    # Get data out from mpileup analysis
    modified_sample_data = data_retriever(m_full_analysis)
    unmodified_sample_data = data_retriever(u_full_analysis)
    shape_data = {'positions':modified_sample_data['positions']}

    # Make overlay plots
    overlay_plot_generator(modified_sample_data ,unmodified_sample_data, output_directory)

    # Subtract background mutation rate
    bg_corr_mut_rate = []
    for i in range(0, len(modified_sample_data["mut_rates"])):
        bg_corr_mut_rate.append(modified_sample_data["mut_rates"][i] - unmodified_sample_data["mut_rates"][i])
    # Add to shape_data dict
    shape_data["reactivity"] = bg_corr_mut_rate

    # Determine uncertainties in mut_rate (stderr)
    shape_stderr = []
    for i in range(0, len(modified_sample_data["mut_rates"])):
        try:
            # Determine independent uncertainty
            mod_stderr = math.sqrt(modified_sample_data["mut_rates"][i]) / math.sqrt(modified_sample_data["depth"][i])
            umod_stderr = math.sqrt(unmodified_sample_data["mut_rates"][i]) / math.sqrt(unmodified_sample_data["depth"][i])
            # Determine propagated uncertainty
            shape_stderr.append(math.sqrt(mod_stderr**2 + umod_stderr**2))
        # If no sequencing depth, no confidence in results at x position
        except ZeroDivisionError:
            mod_stderr = 100
            umod_stderr = 100
            shape_stderr.append(math.sqrt(mod_stderr**2 + umod_stderr**2))
    # Add to shape_data dict
    shape_data["shape_err"] = shape_stderr

    # Make a reactivity plot
    plt.errorbar(
        shape_data["positions"], 
        shape_data["reactivity"], 
        yerr=shape_data["shape_err"], 
        fmt='-', 
        ecolor='black',
        color="orange", 
        capsize=4
        )
    plt.title('Reactivity Plot')
    plt.xlabel('Position')
    plt.ylabel('bg corrected mutation rate + stderr\n(mutations per read)')
    plt.xlim(0, max(shape_data["positions"]))
    plt.ylim(0, max(shape_data["reactivity"])+0.2)
    plt.tight_layout()
    plt.savefig(f"{output_directory}/Reactivity_plot.png")
    plt.clf()

    return shape_data