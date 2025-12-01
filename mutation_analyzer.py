from Bio.Seq import Seq

def clean_sequence(seq):
    seq = seq.upper().replace(" ", "").replace("\n", "")
    valid_bases = {"A", "T", "C", "G"}
    cleaned = "".join([base for base in seq if base in valid_bases])
    return cleaned

# --- User Input ---
original = input("Enter original DNA sequence: ")
mutated = input("Enter mutated DNA sequence: ")

original = clean_sequence(original)
mutated = clean_sequence(mutated)

print("\n--- CLEANED SEQUENCES ---")
print("Original:", original)
print("Mutated:", mutated)
def find_mutation(original, mutated):
    # If lengths are same → substitution
    if len(original) == len(mutated):
        for i in range(len(original)):
            if original[i] != mutated[i]:
                return ("substitution", i, original[i], mutated[i])
        return ("no mutation", None, None, None)

    # If mutated is longer → insertion
    elif len(mutated) > len(original):
        for i in range(len(original)):
            if original[i] != mutated[i]:
                return ("insertion", i, "-", mutated[i])
        return ("insertion", len(original), "-", mutated[len(original)])

    # If mutated is shorter → deletion
    else:  
        for i in range(len(mutated)):
            if original[i] != mutated[i]:
                return ("deletion", i, original[i], "-")
        return ("deletion", len(mutated), original[len(mutated)], "-")

mutation_type, position, old_base, new_base = find_mutation(original, mutated)

print("\n--- MUTATION DETECTED ---")
print("Type:", mutation_type)
print("Position:", position)
print("Original base:", old_base)
print("Mutated base:", new_base)

def translate_dna_to_protein(dna):
    # Trim DNA length so it's divisible by 3 (avoid errors)
    trim_len = len(dna) - (len(dna) % 3)
    dna = dna[:trim_len]
    protein = str(Seq(dna).translate())
    return protein

protein_original = translate_dna_to_protein(original)
protein_mutated = translate_dna_to_protein(mutated)

print("\n--- PROTEIN TRANSLATION ---")
print("Original Protein:", protein_original)
print("Mutated Protein:", protein_mutated)

def compare_proteins(prot1, prot2):
    min_len = min(len(prot1), len(prot2))
    for i in range(min_len):
        if prot1[i] != prot2[i]:
            return i, prot1[i], prot2[i]
    if len(prot1) != len(prot2):
        return min_len, prot1[min_len:], prot2[min_len:]
    return None, None, None

aa_position, aa_original, aa_new = compare_proteins(protein_original, protein_mutated)

print("\n--- AMINO ACID CHANGE ---")
print("Position:", aa_position)
print("Original AA:", aa_original)
print("Mutated AA:", aa_new)

def classify_mutation(mutation_type, aa_original, aa_new):
    if mutation_type == "no mutation":
        return "No impact"

    # Stop codon mutation
    if aa_new == "*" or aa_original == "*":
        return "Nonsense mutation (protein ends early)"

    # No amino acid change
    if aa_original == aa_new:
        return "Silent mutation (no protein change)"

    # Amino acid changed
    if aa_original != aa_new and aa_original is not None:
        return "Missense mutation (protein structure/function may change)"

    # Insertions or deletions that change length
    return "Frameshift mutation (protein drastically altered)"

impact = classify_mutation(mutation_type, aa_original, aa_new)

print("\n--- MUTATION IMPACT ---")
print("Impact:", impact)

import os

# Create outputs folder if not exists
if not os.path.exists("outputs"):
    os.makedirs("outputs")

report_path = "outputs/mutation_report.txt"

with open(report_path, "w", encoding="utf-8") as report:
    report.write("=== DNA MUTATION ANALYSIS REPORT ===\n\n")
    report.write(f"Original DNA: {original}\n")
    report.write(f"Mutated DNA:  {mutated}\n\n")
    report.write("--- Mutation Details ---\n")
    report.write(f"Type: {mutation_type}\n")
    report.write(f"DNA Position: {position}\n")
    report.write(f"Base Change: {old_base} → {new_base}\n\n")
    report.write("--- Protein Translation ---\n")
    report.write(f"Original Protein: {protein_original}\n")
    report.write(f"Mutated Protein:  {protein_mutated}\n\n")
    report.write("--- Amino Acid Change ---\n")
    report.write(f"Position: {aa_position}\n")
    report.write(f"AA Change: {aa_original} → {aa_new}\n\n")
    report.write("--- Biological Impact ---\n")
    report.write(f"{impact}\n\n")
    report.write("=== End of Report ===\n")

print("\nReport generated successfully!")
print(f"File saved at: {report_path}")
