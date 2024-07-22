from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import csv

def calculate_properties(sequence):
    # Remove any gaps (dashes) in the sequence
    cleaned_sequence = sequence.replace('-', '')

    # Analyze the protein sequence, assuming cleaned_sequence is not empty
    if cleaned_sequence:
        analysis = ProteinAnalysis(cleaned_sequence)
        
        # Calculate properties
        protein_length = len(cleaned_sequence)
        hydrophobicity = analysis.gravy()
        charge_at_pH7 = analysis.charge_at_pH(7.0)
        isoelectric_point = analysis.isoelectric_point()
        
        return protein_length, hydrophobicity, charge_at_pH7, isoelectric_point
    else:
        return None

def process_fasta(fasta_file, output_file):
    with open(fasta_file, 'r') as file, open(output_file, 'w', newline='') as out_file:
        writer = csv.writer(out_file, delimiter='\t')
        # Write headers
        writer.writerow(['id', 'length', 'hydrophobicity', 'charge', 'isoelectric_point'])
        
        # Read the fasta file
        for record in SeqIO.parse(file, 'fasta'):
            properties = calculate_properties(str(record.seq))
            if properties:
                # Write the calculated properties to the output file
                writer.writerow([record.id] + list(properties))
            else:
                # Handle sequences that are empty or contain only invalid characters
                writer.writerow([record.id, 'N/A', 'N/A', 'N/A', 'N/A'])

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python script.py input.fasta output.tsv")
    else:
        fasta_file = sys.argv[1]
        output_file = sys.argv[2]
        process_fasta(fasta_file, output_file)
