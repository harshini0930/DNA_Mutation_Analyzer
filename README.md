# DNA_Mutation_Analyzer
A beginner-friendly Python tool to detect DNA mutations, translate sequences into proteins, and predict the biological impact of the mutations.

## Features
Detects Substitution, Insertion, Deletion mutations
Translates DNA sequences to protein sequences
Identifies amino acid changes
Predicts mutation impact: Silent, Missense, Nonsense, Frameshift
Generates a professional report (`mutation_report.txt`)

## Requirements
 Python 3.x
BioPython (`pip install biopython`)

## Usage
Clone or download this repository.
Open terminal in the project folder.
Run the script:
    python mutation_analyzer.py
Enter the original DNA sequence when prompted.
Enter the mutated DNA sequence when prompted.
Check the outputs/mutation_report.txt file for the full report.

## Example
Input:
    Original DNA: ATGGAATCT
    Mutated DNA:  ATGGGATCT

Output:
    === DNA MUTATION ANALYSIS REPORT ===

Original DNA: ATGGAATCT
Mutated DNA:  ATGGGATCT

--- Mutation Details ---
Type: substitution
DNA Position: 4
Base Change: A → G

--- Protein Translation ---
Original Protein: MES
Mutated Protein:  MGS

--- Amino Acid Change ---
Position: 1
AA Change: E → G

--- Biological Impact ---
Missense mutation (protein structure/function may change)

=== End of Report ===

## Future Improvements
Add a GUI with Tkinter for easy input.
Add visualizations for amino acid changes.
Integrate real mutation databases (ClinVar, dbSNP).
Add ML prediction for pathogenicity of mutations.
Integrate real mutation databases (ClinVar, dbSNP).
Add ML prediction for pathogenicity of mutations.
