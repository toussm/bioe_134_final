import traceback
import re
import csv
import time
from dnachisel import *
from transcriptDesigner.transcript_designer import TranscriptDesigner

def parse_fasta(fasta_file):
    """
    Parses the FASTA file to extract gene names and protein sequences.
    """
    sequences = {}
    current_gene = None
    current_sequence = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_gene:
                    sequences[current_gene] = ''.join(current_sequence)
                gene_name = None
                parts = line.split()
                for part in parts:
                    if part.startswith("GN="):
                        gene_name = part.split("=")[1]
                        break
                if not gene_name:
                    gene_name = line.split('|')[2].split(' ')[0]
                current_gene = gene_name
                current_sequence = []
            else:
                current_sequence.append(line)
        if current_gene:
            sequences[current_gene] = ''.join(current_sequence)
    
    return sequences

def benchmark_proteome(fasta_file):
    """
    Benchmarks the proteome using TranscriptDesigner.
    """
    designer = TranscriptDesigner()
    designer.initiate()

    proteome = parse_fasta(fasta_file)
    successful_results = []
    error_results = []
    i = 0
    
    for gene, protein in proteome.items():
        try:
            print(f"Processing gene: {gene} with protein sequence: {protein[:30]}...")
            transcript = designer.run(protein)
            successful_results.append({
                'gene': gene,
                'protein': protein,
                'transcript': transcript[0],
                'summary' : transcript[1]
            })
            # if i == 2:
            #     break
        except Exception as e:
            error_results.append({
                'gene': gene,
                'protein': protein,
                'error': f"Error: {str(e)}\nTraceback: {traceback.format_exc()}"
            })
            if i == 2:
                break

        i+=1
    return successful_results, error_results

def analyze_errors(error_results):
    """
    Write the error analysis to a text file.
    """
    error_summary = {}
    with open('error_summary.txt', 'w') as f:
        for error in error_results:
            error_message = error['error'].split("\n")[0]
            error_summary[error_message] = error_summary.get(error_message, 0) + 1
            f.write(f"Gene: {error['gene']}\n{error['error']}\n\n")
    
    return error_summary

def parse_constraints_summary(summary, gene, protein, transcript_dna, validation_failures):

    constraint_aggregates = {
        "Forbidden Patterns": [],
        "Hairpins": [],
        "GC Content": [],
        "Rare Codons": []
    }

    # Check for failed AvoidPattern constraints
    avoid_pattern_failures = re.findall(
        r'FAIL ┍ AvoidPattern\[.*?\]\(pattern:(\S+)\).*?Pattern found at positions (.*?)\n',
        summary,
        re.DOTALL
    )
    # if len(avoid_pattern_failures) > len(transcript_dna) // 300:
    #     constraint_aggregates["Forbidden Patterns"].append("True")
    for pattern, locations in avoid_pattern_failures:
        formatted_locations = locations.replace('\n', ' ').replace('"', "'").strip()
        constraint_aggregates["Forbidden Patterns"].append(f"{pattern} at {formatted_locations}")
    
    
    avoid_rare_codons = re.findall(
        r'FAIL ┍ AvoidRareCodons\[.*?\]\n.*?Pattern found at positions (.*?)\n',
        summary,
        re.DOTALL
    )

    for pattern, locations in avoid_rare_codons:
        formatted_locations = locations.replace('\n', ' ').replace('"', "'").strip()
        constraint_aggregates["Rare Codons"].append(f"{formatted_locations}")

    # Check for failed AvoidHairpins constraints
    failed_hairpins = re.findall(
        r'FAIL ┍ AvoidHairpins\[.*?\].*?Locations: (.*?)\n',
        summary,
        re.DOTALL
    )
    for hairpin_locations in failed_hairpins:
        formatted_hairpin = hairpin_locations.replace('\n', ' ').replace('"', "'").strip()
        constraint_aggregates["Hairpins"].append(formatted_hairpin)

    # Check for failed EnforceGCContent constraints
    gc_content_failures = re.findall(
        r'FAIL ┍ EnforceGCContent\[.*?\]\(mini:(\S+), maxi:(\S+)\)',
        summary
    )
    for mini, maxi in gc_content_failures:
        constraint_aggregates["GC Content"].append(f"Out of bounds: {mini}-{maxi}")

    # Add aggregated failures to validation_failures
    for constraint_type, failures in constraint_aggregates.items():
        if failures:
            aggregated_sites = "; ".join(failures)
            validation_failures.append({
                'gene': gene,
                'protein': protein,
                'cds': transcript_dna,
                'site': f"{constraint_type} detected: {aggregated_sites}"
            })

    return validation_failures

def validate_transcripts(successful_results):
    """
    Validate the successful transcripts using various checkers, now including CodonChecker.
    """

    validation_failures = []
    for result in successful_results:
        cds = ''.join(result['transcript'])
        try:
            # Check if CDS length is a multiple of 3
            if len(cds) % 3 != 0:
                raise ValueError("CDS length is not a multiple of 3.")

            # Verify that the translated protein matches the original protein
            original_protein = result['protein']
            translated_protein = translate(cds[:-1])
            mismatches = []
            min_length = min(len(original_protein), len(translated_protein))

            for i in range(min_length):
                if original_protein[i] != translated_protein[i]:
                    mismatches.append({
                        'position': i + 1,  # Residue positions are 1-based
                        'original': original_protein[i],
                        'translated': translated_protein[i]
                    })
            # print(mismatches)
            if original_protein != translated_protein:
                raise ValueError(f"Translation mismatch: Original {original_protein}, Translated {translated_protein}")

            # Ensure CDS starts with valid start codon and ends with stop codon
            if not (cds.startswith(("ATG", "GTG", "TTG", "CTG")) and cds.endswith(("TAA", "TGA", "TAG"))):
                raise ValueError("CDS does not start with a valid start codon or end with a valid stop codon.")
        except ValueError as e:
            validation_failures.append({
                'gene': result['gene'],
                'protein': result['protein'],
                'cds': cds,
                'site': f"Translation or completeness error: {str(e)}"
            })
            continue

        # Validate against hairpins, forbidden sequences, and codons
        parse_constraints_summary(validation_failures=validation_failures, summary=result['summary'], gene=result['gene'], protein=result['protein'], transcript_dna=result['transcript'])
    return validation_failures

def write_validation_report(validation_failures):
    """
    Writes validation results to a TSV file.
    """
    with open('validation_failures.tsv', 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['gene', 'protein', 'cds', 'site'])
        for failure in validation_failures:
            writer.writerow([failure['gene'], failure['protein'], failure['cds'], failure['site']])

def generate_summary(total_genes, parsing_time, execution_time, errors_summary, validation_failures):
    """
    Generates a streamlined summary report categorizing validation failures by checker.
    """
    total_validation_failures = len(validation_failures)
    
   # Categorize failures by checker type
    checker_failures = {
        'Forbidden Sequence Checker': 0,
        'Hairpin Checker': 0,
        'GC Content Checker': 0,
        'Rare Codons': 0,
        'Translation/Completeness Checker': 0
    }

    # Increment the appropriate checker category based on the failure site
    for failure in validation_failures:
        site = failure['site']
        if "Forbidden Patterns" in site:
            checker_failures['Forbidden Sequence Checker'] += 1
        elif "Hairpins" in site:
            checker_failures['Hairpin Checker'] += 1
        elif "GC Content" in site:
            checker_failures['GC Content Checker'] += 1
        elif "Rare Codons" in site:
            checker_failures['Rare Codons'] += 1
        elif "Translation or completeness error" in site:
            checker_failures['Translation/Completeness Checker'] += 1

    # Generate the summary report
    with open('summary_report.txt', 'w') as f:
        f.write(f"Total genes processed: {total_genes}\n")
        f.write(f"Parsing runtime: {parsing_time:.2f} seconds\n")
        f.write(f"Execution runtime: {execution_time:.2f} seconds\n")
        f.write(f"Total exceptions: {sum(errors_summary.values())}\n")

        if errors_summary:
            f.write(f"\nTop 3 most common exceptions:\n")
            for error, count in sorted(errors_summary.items(), key=lambda x: x[1], reverse=True)[:3]:
                f.write(f"- {error}: {count} occurrences\n")
        else:
            f.write("No exceptions encountered.\n")

        f.write(f"\nTotal validation failures: {total_validation_failures}\n")

        # Categorize validation failures by checker
        f.write("\nValidation Failures by Checker:\n")
        for checker, count in checker_failures.items():
            f.write(f"- {checker}: {count} occurrences\n")

def run_benchmark(fasta_file):
    """
    Runs the complete benchmark process: parsing, running TranscriptDesigner, validating, and generating reports.
    """
    start_time = time.time()
    
    # Benchmark the proteome
    parsing_start = time.time()
    successful_results, error_results = benchmark_proteome(fasta_file)
    parsing_time = time.time() - parsing_start
    
    # Analyze and log errors
    errors_summary = analyze_errors(error_results)
    
    # Validate the successful transcripts
    validation_start = time.time()
    validation_failures = validate_transcripts(successful_results)
    execution_time = time.time() - validation_start

    # Write validation and error reports
    write_validation_report(validation_failures)
    
    # Generate the summary report
    total_genes = len(successful_results) + len(error_results)
    generate_summary(total_genes, parsing_time, execution_time, errors_summary, validation_failures)

if __name__ == "__main__":
    fasta_file = "benchmarker/UP000002311_559292.fasta"
    run_benchmark(fasta_file)
