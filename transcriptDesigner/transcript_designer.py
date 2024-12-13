from dnachisel import *
from transcriptDesigner.one_to_one_revTranslate import amino_acid_to_codon

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

    
    def optimize_transcript(self, dna_sequence: str):
        # Define the optimization problem
        problem = DnaOptimizationProblem(
            sequence=dna_sequence,
            constraints=[
                EnforceGCContent(mini=0.2, maxi=0.4),
                AvoidHairpins(stem_size=10),  # Avoid hairpins in the sequence
                AvoidPattern("GAATTC"),  # Avoid EcoRI restriction site
                AvoidPattern("GGATCC"),  # Avoid BamHI restriction site
                AvoidPattern("AAAAAAAAA"),  # Avoid long homopolymer stretches
                AvoidPattern("TATAAA"),  # Avoid TATA box
            ],
            objectives=[
                CodonOptimize(species='s_cerevisiae'),  # Optimize for S. cerevisiae codon usage
            ]
        )

        # Solve the problem
        problem.optimize()
        print(problem.sequence)

        # Ensure the solution is valid
        if problem.all_constraints_pass():
            return problem.sequence  # Optimized DNA sequence
        else:
            print(problem.constraints_text_summary())
            raise ValueError("Optimization failed: constraints could not be satisfied.")
        
    def run(self, protein: str):

        #Check for stop codon. Add if not there
        if protein[-1] != '*':
            protein += '*'
        
        simple_dna_seq = self.simpleReverseTranslate(protein)
        optimized_transcript = self.optimize_transcript(simple_dna_seq)
        
        return optimized_transcript
    
if __name__ == "__main__":
    # Example usage of TranscriptDesigner
    peptide = "MGLILRWKEKKQLSSKQNAQKSRKPANTSFRQQRLKAWQPILSPQSVLPL"
    
    designer = TranscriptDesigner()
    designer.initiate()

    transcript = designer.run(peptide)
    translation = translate(transcript)

    print(transcript)
    print(translation)
       