#!/usr/bin/env python3

import sys
import pandas as pd
import matplotlib.pyplot as plt

def generate_coverage_plot(input_file, output_file):
    """
    Generate a coverage plot from a coverage file.

    Args:
        input_file (str): Path to the input coverage file (e.g., `${barcode}_coverage.txt`).
        output_file (str): Path to save the output plot (e.g., `${barcode}_coverage_plot.png`).
    """
    # Read the coverage file into a pandas DataFrame
    try:
        # Assuming the file is a TSV (tab-separated values) file
        df = pd.read_csv(input_file, sep='\t')
    except Exception as e:
        print(f"Error reading input file: {e}")
        sys.exit(1)

    # Check if the required columns are present
    # Based on `samtools coverage` output format
    required_columns = ['#rname', 'startpos', 'endpos', 'coverage']
    if not all(col in df.columns for col in required_columns):
        print(f"Input file must contain the following columns: {required_columns}")
        sys.exit(1)

    # Calculate the midpoint of each region for plotting
    df['midpos'] = (df['startpos'] + df['endpos']) / 2

    # Plot the coverage
    plt.figure(figsize=(10, 6))
    plt.plot(df['midpos'], df['coverage'], label='Coverage', color='blue', marker='o', linestyle='-')
    plt.title('Genome Coverage Plot')
    plt.xlabel('Genomic Position')
    plt.ylabel('Coverage Depth')
    plt.grid(True)
    plt.legend()

    # Save the plot to a file
    plt.savefig(output_file)
    print(f"Coverage plot saved to {output_file}")

if __name__ == "__main__":
    # Parse command-line arguments
    if len(sys.argv) != 3:
        print("Usage: python3 generate_coverage_plot.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Generate the coverage plot
    generate_coverage_plot(input_file, output_file)