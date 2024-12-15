from dnachisel import *
from transcriptDesigner.one_to_one_revTranslate import amino_acid_to_codon
import warnings
from Bio import BiopythonWarning

warnings.simplefilter('ignore', BiopythonWarning)

class TranscriptDesigner:
    def __init__(self):
        self.aminoAcidToCodon = {}

    def initiate(self) -> None:
        """
        Initializes the codon table, the RBS chooser and the checkers.
        """
        self.aminoAcidToCodon = amino_acid_to_codon('yeast_codon_usage_table.tsv')
    
    def simpleReverseTranslate(self, peptide: str):
        dna_seq = ""
        for aa in peptide:
            codon = self.aminoAcidToCodon[aa]
            dna_seq += codon
        return dna_seq

    
    def optimize_transcript(self, dna_sequence: str, peptide: str):
        # Define the optimization problem
        problem = DnaOptimizationProblem(
            sequence=dna_sequence,
            constraints=[
                EnforceGCContent(mini=0.2, maxi=0.6),
                AvoidHairpins(stem_size=10, hairpin_window=30),  # Avoid hairpins in the sequence
                AvoidPattern("GAATTC"),  # Avoid EcoRI restriction site
                AvoidPattern("GGATCC"),  # Avoid BamHI restriction site
                AvoidPattern("AAAAAAAAA"),  # Avoid long homopolymer stretches
                AvoidPattern("TATAAA"),  # Avoid TATA box
                AvoidPattern("AAGCTT"),  # Avoid HindIII restriction site
                AvoidPattern("GCGGCCGC"), # Avoid NotI restriction site
                AvoidPattern("GTAAGT"), # Avoid intron splice sites
                AvoidPattern("GTGAGT"), # Avoid intron splice sites
                AvoidRareCodons(species="s_cerevisiae", min_frequency=0.2),
                EnforceTranslation(translation=peptide)

            ],
            objectives=[
                EnforceGCContent(mini=0.2, maxi=0.6),
                AvoidPattern("GAATTC"),  # Avoid EcoRI restriction site
                AvoidPattern("GGATCC"),  # Avoid BamHI restriction site
                AvoidPattern("AAAAAAAAA"),  # Avoid long homopolymer stretches
                AvoidPattern("TATAAA"),  # Avoid TATA box
                AvoidPattern("AAGCTT"),  # Avoid HindIII restriction site
                AvoidHairpins(stem_size=10, hairpin_window=30),  # Avoid hairpins in the sequence
                AvoidPattern("GCGGCCGC"), # Avoid NotI restriction site
                AvoidPattern("GTAAGT"), # Avoid intron splice sites
                AvoidPattern("GTGAGT"), # Avoid intron splice sites
                AvoidRareCodons(species="s_cerevisiae", min_frequency=0.2),
                EnforceTranslation(translation=peptide),
                CodonOptimize(species='s_cerevisiae')  # Optimize for S. cerevisiae codon usage
            ]
        )

        # Solve the problem
        problem.resolve_constraints()
        problem.optimize()
        return problem
        
    def run(self, protein: str):

        protein += '*'
        
        simple_dna_seq = self.simpleReverseTranslate(protein)
        optimized_transcript = self.optimize_transcript(simple_dna_seq, protein)
        
        return optimized_transcript.sequence, optimized_transcript.constraints_text_summary()
    
if __name__ == "__main__":
    # Example usage of TranscriptDesigner
    peptide = "MGLILRWKEKKQLSSKQNAQKSRKPANTSFRQQRLKAWQPILSPQSVLPL"
    
    designer = TranscriptDesigner()
    designer.initiate()

    transcript = designer.run(peptide)[0]
    translation = translate(transcript)

    print(transcript)
    print(translation)
       