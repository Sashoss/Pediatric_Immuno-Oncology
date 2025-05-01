
import os
import numpy as np

def read_rmsf_file(file_path):
    """
    Reads an RMSF .xvg file and extracts the numeric data.
    Assumes RMSF values are in the second column and first column is atom numbers.
    """
    atom_numbers = []
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            if not line.startswith(('@', '#')):  
                parts = line.strip().split()
                atom_numbers.append(int(parts[0]))  
                data.append(float(parts[1]))  
    return atom_numbers, np.array(data)

def calculate_position_wise_quadratic_mean(runs_folder, runs=['Run_MD_chain1', 'Run_MD_chain2', 'Run_MD_chain3'], protein_chain="SPP1"):
    """
    Calculates the quadratic mean of RMSF values at each position across multiple runs.
    Assumes atom numbers are consistent across all runs.
    """
    all_rmsf = []
    atom_numbers = None
    sample=runs_folder.split("/")[2]
    for run in runs:
        chain_id = run.split("_")[2]
        rmsf_file = os.path.join(runs_folder, run, 'rmsf_'+sample+"_"+chain_id+"_"+protein_chain+'.xvg')
        if os.path.exists(rmsf_file):
            atoms, rmsf_data = read_rmsf_file(rmsf_file)
            if atom_numbers is None:
                atom_numbers = atoms
            elif atom_numbers != atoms:
                raise ValueError("Atom numbers do not match across runs")
            all_rmsf.append(rmsf_data)
        else:
            print(f"File not found: {rmsf_file}")
            return None, None

    all_rmsf = np.vstack(all_rmsf)
    quadratic_mean_per_position = np.sqrt(np.mean(all_rmsf**2, axis=0))
    return atom_numbers, quadratic_mean_per_position

def write_xvg(file_path, atom_numbers, data, title, x_label="Atom Number", y_label="Quadratic Mean RMSF"):
    """
    Writes data to an .xvg file with appropriate headers.
    """
    with open(file_path, 'w') as f:
        f.write(f"@    title \"{title}\"\n")
        f.write(f"@    xaxis  label \"{x_label}\"\n")
        f.write(f"@    yaxis  label \"{y_label}\"\n")
        f.write(f"@TYPE xy\n")
        for atom, value in zip(atom_numbers, data):
            f.write(f"{atom} {value:.6f}\n")

def main():
    base_dir = "." 
    folders = ['out/23C3_v1', 'out/Hu23C3_v1']
    results = {}

    for protein_chain in ["SPP1", "Heavy", "Light"]:
        for folder in folders:
            runs_folder = os.path.join(base_dir, folder, 'md')
            atom_numbers, quadratic_mean = calculate_position_wise_quadratic_mean(runs_folder, protein_chain=protein_chain)
            if quadratic_mean is not None:
                results[folder] = (atom_numbers, quadratic_mean)
                output_file = os.path.join(base_dir, f"{folder}_{protein_chain}_quadratic_mean.xvg")
                write_xvg(output_file, atom_numbers, quadratic_mean, title=f"Quadratic Mean RMSF for {folder}")
                print(f"Written quadratic mean to {output_file}")

if __name__ == "__main__":
    main()

