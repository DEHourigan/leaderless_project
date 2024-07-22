import os
import pandas as pd
import argparse

def append_assembly_data(protein_properties_file, bakta_out_dir, output_file):
    # Define the new column names based on the provided schema
    column_names = ['sequence_id', 'type', 'start', 'stop', 'strand', 'id', 'gene', 'product', 'dbxrefs']
    
    # Load the protein properties file, ensuring 'id' is treated as a string
    protein_properties = pd.read_csv(protein_properties_file, sep='\t', dtype={'id': str})
    
    # Iterate over all tsv files in the bakta_out directory
    for filename in os.listdir(bakta_out_dir):
        if filename.endswith('.tsv') and not filename.endswith('hypotheticals.tsv'):
            file_path = os.path.join(bakta_out_dir, filename)
            try:
                # Load the assembly data, skipping lines that start with '#' and without headers
                assembly_data = pd.read_csv(file_path, sep='\t', comment='#', header=None, names=column_names)
                
                # Ensure 'id' is in the dataframe
                if 'id' in assembly_data.columns:
                    # Merge the data on 'id', adding the assembly information
                    # Use 'how=left' to keep all entries from protein_properties and only add matching from assembly_data
                    protein_properties = pd.merge(protein_properties, assembly_data, on='id', how='left')
            except pd.errors.EmptyDataError:
                print(f"No data in {file_path}, skipping.")
            except Exception as e:
                print(f"Error processing {file_path}: {e}")

    # Save the updated protein properties file
    protein_properties.to_csv(output_file, sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Append assembly data to protein properties file.')
    parser.add_argument('protein_properties_file', type=str, help='Path to the protein properties TSV file.')
    parser.add_argument('bakta_out_dir', type=str, help='Directory path containing .tsv files.')
    parser.add_argument('output_file', type=str, help='Path to the output file where updated TSV will be saved.')
    args = parser.parse_args()

    append_assembly_data(args.protein_properties_file, args.bakta_out_dir, args.output_file)
