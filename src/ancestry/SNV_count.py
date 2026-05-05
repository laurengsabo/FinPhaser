import pandas as pd
import numpy as np
import glob
import os
import sys

def process_ancestry_files(target_dir):
    # 1. Setup paths
    search_path = os.path.join(target_dir, "*.posterior")
    files = glob.glob(search_path)
    
    # Create the new output directory
    output_dir = os.path.join(os.getcwd(), "ancestry_summary_report")
    os.makedirs(output_dir, exist_ok=True)

    if not files:
        print(f"Error: No .posterior files found in directory: {target_dir}")
        return

    summary_data_list = []

    print(f"Found {len(files)} files. Starting processing...")

    for file_path in files:
        file_name = os.path.basename(file_path)
        
        # Load the data
        df = pd.read_csv(file_path, sep='\t')
        
        # Identify the probability columns (skipping chrom and position)
        prob_cols = df.columns[2:]
        
        # Calculate the most likely state for each row
        best_indices = np.argmax(df[prob_cols].values, axis=1)
        df['called_state'] = prob_cols[best_indices]
        
        # Record the counts for the master summary
        counts = df['called_state'].value_counts()
        counts.name = file_name
        summary_data_list.append(counts)
        
        # Define output path within the new directory
        output_csv = os.path.join(output_dir, f"calls_{file_name}.csv")
        
        # Write individual file with summary information in the header
        with open(output_csv, 'w') as f:
            f.write(f"# File: {file_name}\n")
            f.write("# Summary of state calls for this individual:\n")
            for state, count in counts.items():
                f.write(f"#   {state}: {count}\n")
            f.write("#" + "-"*40 + "\n")
            
            # Write the data
            df[['chrom', 'position', 'called_state']].to_csv(f, index=False)
        
        print(f"  Finished: {file_name}")

    # Create the Master Summary Report
    summary_df = pd.concat(summary_data_list, axis=1).T.fillna(0).astype(int)
    
    # Calculate global totals across all files
    global_totals = summary_df.sum(axis=0).sort_values(ascending=False)
    
    summary_report_path = os.path.join(output_dir, "master_ancestry_summary.csv")
    
    with open(summary_report_path, 'w') as f:
        f.write("# GLOBAL ANCESTRY SUMMARY (Aggregated over all files)\n")
        f.write("# Total SNVs called per state:\n")
        for state, total in global_totals.items():
            f.write(f"#   {state}: {total}\n")
        f.write("#" + "-"*40 + "\n")
        
        # Write the per-file summary table
        summary_df.to_csv(f)

    print("\nProcessing complete!")
    print(f"All outputs saved to the directory: {output_dir}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        process_ancestry_files(sys.argv[1])
    else:
        print("Usage: python script_name.py <directory_with_posterior_files>")