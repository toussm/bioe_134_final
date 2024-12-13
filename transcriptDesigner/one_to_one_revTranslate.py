import csv

def amino_acid_to_codon(file_path):
    codon_dict = {}

    # Open and read the TSV file
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')

        for row in reader:
            codon = row['Codon']
            amino_acid = row['Amino Acid']
            frequency = float(row['Frequency'])

            # Skip rows with missing codons or amino acids
            if not codon or not amino_acid:
                continue

            # Check if this codon has the highest frequency for the amino acid
            if amino_acid not in codon_dict or frequency > codon_dict[amino_acid]['Frequency']:
                codon_dict[amino_acid] = {'Codon': codon, 'Frequency': frequency}

    # Convert the dictionary to a simple mapping of amino acid to codon
    amino_acid_to_codon = {aa: data['Codon'] for aa, data in codon_dict.items()}
    return amino_acid_to_codon

# Example usage:
file_path = "yeast_codon_usage_table.tsv"  # Replace with your file path
result = amino_acid_to_codon(file_path)