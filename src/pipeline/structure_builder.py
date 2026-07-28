import argparse, time, sys, subprocess
import numpy as np
from pathlib import Path

class UnionFind:
    # Disjoint Set data structure
    def __init__(self, size):
        # Each element is its own representative
        self.parent = list(range(size))

    def find(self, i):
        # Is 'i' a root/representative of a group?
        if self.parent[i] == i:
            return i
        # Else, recursively find root/representative
        else:
            return self.find(self.parent[i])
        
    def unite(self, i, j):
        # Make a connection between two representatives
        # Representatives i and j
        irep = self.find(i)
        jrep = self.find(j)

        # The representative of i's set is now the representative of j's set
        self.parent[irep] = jrep


def predict_structs(sequence, constraints=None):
    # TODO
    pass

def read_bracket_dot(dot_bracket):
    data = {}
    with open(dot_bracket, "r") as f:
        lines = []
        for line in f:
            if line.strip():
                lines.append(line.strip())
        # Process lines
        i = 0
        count = 1
        while i < len(lines):
            # Debug:
            # print(f"{i} / {len(lines)}")
            # print(count)
            # Entries begin with ">"
            if str(lines[i]).startswith(">"):
                # capture header - contains structure ∆G in kcal/mol
                header = lines[i]
                # First structure also has sequence
                if i == 0:
                    sequence = lines[i+1]
                    structure = lines[i+2]
                    # Add data to output
                    data["sequence"] = sequence
                    data[f"structure_{count}"] = {"header":header, "structure":structure}
                    # Increase iteration and captured sequence count
                    count += 1
                    i += 3
                # Other structures do not include the sequence line
                else:
                    structure = lines[i+1]
                    # Add data to output
                    data[f"structure_{count}"] = {"header":header, "structure":structure}
                    # Increase iteration and captured sequence count
                    count += 1
                    i += 2
    # Return data dict
    return data

def db_to_bin(db_data):
    out_str = ""
    # Convert dot bracket to binary encoding of base pair True / False
    for position in range(0,len(db_data)):
        # Open or Close base pair is represented as 1
        if db_data[position] == "(":
            out_str += "1"
        elif db_data[position] == ")":
            out_str += "1"
        # Unpaired base is represented as 0
        elif db_data[position] == ".":
            out_str += "0"
        # Should probably handle this better
        else:
            print(f"[!] Unknown Symbol! - '{db_data[position]}'")
    return out_str

def compare_structures(dot_bracket_data:dict):
    structure_bins = {}
    # Get all dot bracket structures + convert to binary encoding for comparison
    for key in dot_bracket_data.keys():
        if key == "sequence":
            continue
        else:
            structure_bins[key] = db_to_bin(dot_bracket_data[key]["structure"])

    # Pairwise score dict
    similarity_matrix = {}

    # Set one structure as reference
    for ref_key, ref_bin in structure_bins.items():
        # Initial matrix entry for ref_bin
        similarity_matrix[ref_key] = {}

        # Compare against all other structures
        for alt_key, alt_bin in structure_bins.items():
            if ref_key == alt_key:
                # Skip self-comparison
                continue

            # Initial comparison score
            score = 0
            for pos in range(len(ref_bin)):
                # +1 for matches
                if ref_bin[pos] == alt_bin[pos]:
                    score += 1
            similarity_matrix[ref_key][alt_key] = score

    # Find and print most similar pairs:
    pairs = {}
    for ref_key, comparisons in similarity_matrix.items():
        if not comparisons:
            continue
        # Get highest score
        best_match = max(comparisons, key=comparisons.get)
        best_match_score = comparisons[best_match]
        all_scores = []
        for struct, scores in comparisons.items():
            all_scores.append((struct, scores))
        # Arbitrary 80% match threshold
        if round(best_match_score / len(ref_bin),2) >= 0.8:
            print(f"[*] '{ref_key}' is most similar to '{best_match}'\n\tScore: {best_match_score} / {len(ref_bin)}; {best_match_score / len(ref_bin) * 100:1.0f}% >> {all_scores}\n")
            # Collapse down to 0-based indexing and retain best match
            pairs[int(ref_key[-1])-1] = int(best_match[-1])-1
        else:
            # If best match does not meet threshold, no relation to other structures (best match is itself)
            print(f"[*] '{ref_key}' is most similar to '{ref_key}'\n\t{all_scores}\n")
            pairs[int(ref_key[-1])-1] = int(ref_key[-1])-1

    num_structs = len(pairs.keys())
    # Form Disjoint sets / Union Find
    uf = UnionFind(num_structs)
    # Pair up all relevant structures (i.e. max comparison score)
    for mate_1, mate_2 in pairs.items():
        uf.unite(mate_1, mate_2)

    subgroups_master = []
    for index in range(num_structs):
        # What subgroup does structure belong to?
        subgroup = uf.find(index)
        # Put true index in master
        subgroups_master.append(subgroup)
    # Get final .parent list for structure subgroups
    uf.parent = subgroups_master
    # Debug:
    # print(uf.parent)
    # Output a list of disjoint sets with the comparable structures for consensus analysis
    print(f"[-] Number of structure subgroups: {len(set(uf.parent))}\n")
    return uf

def parse_structure(structure):
    pairs = []
    stack = []
    for index, char in enumerate(structure):
        if char == "(":
            stack.append(index)
        elif char == ")":
            if stack: # If stack is empty, evals False
                start = stack.pop()
                pairs.append((start, index))
    return pairs

def consensus_structure(dot_bracket_data:dict, uf):
    structures = {}
    final_structures = {}
    # Init subgrouping
    subgroups = {i:[] for i in set(uf.parent)}

    # Get all dot bracket structures
    for key in dot_bracket_data.keys():
        if key == "sequence":
            continue
        else:
            structures[key] = dot_bracket_data[key]["structure"]

    # For each structure
    for index, struct in enumerate(structures.values()):
        # What subgroup does structure belong to?
        subgroup = uf.find(index)
        # Add structure to subgroup
        subgroups[subgroup].append(struct)
    
    # For each subgroup, get a consensus
    for grouped_structs in subgroups.values():
        # Tracking the 1-based index of subgroup for print messages and output
        subgroup_index = list(subgroups.values()).index(grouped_structs)+1
        # Accumulate majority pairing in structure
        pair_counts = {}
        
        # If there's only one structure in subgroup...
        if len(grouped_structs) == 1:
            consensus = grouped_structs[0]
            print(f"[-] Consensus structure for subgroup {subgroup_index}: {consensus}\n")
            final_structures[f"consensus_struct_{subgroup_index}"] = consensus
            dot_bracket_data['consensus-structures'] = final_structures
            continue

        # Get all symbols at position in each structure
        for struct in grouped_structs:
            for start, end in parse_structure(struct):
                # Sum how many times a specific pair was found
                pair_counts[(start, end)] = pair_counts.get((start, end), 0) + 1
        # Keep pairs wih a frequence of at least 2, then sort
        # x is a key (pair_counts.keys), fetch negative frequency (sorting highest frequencies first), start position is a tiebreaker
        # yields a list of keys
        valid_pairs = sorted([pair for pair, count in pair_counts.items() if count >= 2], key=lambda x: (-pair_counts[x], x[0]))
        print(f"[*] Valid pairs for subgroup {subgroup_index}:\n{valid_pairs}")
        print(f"[*] Frequencies: {[pair_counts[pair] for pair in valid_pairs]}\n")
        # Resolve duplicate pairs
        final_pairs = []
        used_indices = set()

        for i, j in valid_pairs:
            # Prevent multiple pairs from single position
            if i not in used_indices and j not in used_indices:
                # Hence why pairs are sorted by frequency before position
                final_pairs.append((i,j))
                used_indices.add(i)
                used_indices.add(j)   

        # Recapitulate structure
        consensus = list("." * len(grouped_structs[0]))
        for i, j in final_pairs:
            consensus[i] = "("
            consensus[j] = ")"

        print(f"[-] Consensus structure for subgroup {subgroup_index}: {"".join(consensus)}\n")
        final_structures[f"consensus_struct_{subgroup_index}"] = consensus
        dot_bracket_data['consensus-structures'] = final_structures
    return dot_bracket_data

def write_out_consensus(dot_bracket:Path, dot_bracket_data:dict):
        for key, struct in dot_bracket_data['consensus-structures'].items():
            output_file = dot_bracket.parent / f"{key}.dbn"
            rna_name = dot_bracket_data['structure_1']['header'].split(" ")[-1]
            with open(output_file, "w") as f:
                f.write(f">ENERGY = -999 {rna_name}\n{dot_bracket_data['sequence']}\n")
                f.write(f"{"".join(struct)}")

def main():
    # Track evaluation time
    start = time.perf_counter()
    print("[*] Starting structure generator")

    # Parse arguments
    parser = argparse.ArgumentParser(description="Secondary structure prediction (individual and consensus)")
    parser.add_argument("--ref", type=str, help=".fasta file containing the sequence to predict")
    parser.add_argument("--map", type=str, help=".map file to constrain predicted structures from shapemapper2")
    parser.add_argument("--bd", type=str, help=".dbn file generated by RNAStructure - structures to be collapsed to consensus")

    args = parser.parse_args()
    # If no args, print help message
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # TODO - Validate args
    try:
        ref_path = Path(args.ref)
        map_path = Path(args.map)
    except TypeError:
        ref_path = None
        map_path = None
    try:
        bd_path = Path(args.bd)
    except TypeError:
        bd_path = None

    # If --ref and/or --map, generate predicted secondary structures via RNAStructure "Fold" tool
    if ref_path != None and map_path != None:
        # Predict structs with SHAPE constraints
        predict_structs(ref_path, map_path)
    elif ref_path != None and map_path == None:
        predict_structs(ref_path)
    elif ref_path == None and map_path == None:
        pass
    else:
        print("[!] Something went wrong, check the input paths")

    # If --bd, collapse the secondary structures in given bracket-dot to a single consensus structure
    if bd_path != None:
        # Get data
        structure_data = read_bracket_dot(bd_path)
        # Collapse structures to consensus
        compared_structures = compare_structures(structure_data)
        final_data = consensus_structure(structure_data, compared_structures)
        write_out_consensus(bd_path, final_data)
    
    end = time.perf_counter()
    print(f"[-] Structure(s) built in {end - start:.2f} seconds")
if __name__ == "__main__":
    main()