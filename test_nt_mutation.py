"""
Test script for nucleotide-level mutation scanning
"""

from tb_resistance_mvp.pipeline import mutational_scan

# Test sequence (rpoB fragment)
nt_sequence = "ATGGCGTCAACAAAGCAGCTGCTGGCGGTCCATGTGCCGCGTAACACTGGACGAATATCAGTTTGGCCTGAGCTGGGTGAAGATGACGCCGAACATCGCCGAGCGCTGCCTGGTCGATTCGGGCCACACGGTGCAGTTCCCGGCGCGGCTCAAGATGGGCTACGATGAAACCATCGTGAACCAGGCGCTGAGCCGTCCGTGGTTCAAGGGCTGCTACATGGACACCAACAGCCTGGAGATCCGGGCGCAGCCGCTGCACGGGATGTCGTTCTACACGAAGGTGGACGCCTGCATCGGCCTGCCGAACCGGCACACGCAGGTGGAGTGGATGGACGCCTCATTCAAGCGGCTGTACACCCCGGGCCAGAACGTGCTGTCGGAGGCCTGCCGGAAGATCTTCCACGACGGCACCATGGTGCTGCCGCAGGCGAACAGCCGCTGGGACTATGAGAAGGGCACGATGCTGTTCCCGATCAGCCACGTCCAGCGCGACGCGAACTGCCTGGGCACCGAATACTGGATGAAGCCGCGCAGCTTCGCCGTGCTGCAGAACGACGGCAGCATCACCCGCAAGATGCTGGAGCACGTGCCGGTG"

# Translate to protein (manual check)
from tb_resistance_mvp.features import translate_codon

protein = ""
for i in range(0, len(nt_sequence), 3):
    codon = nt_sequence[i:i+3]
    if len(codon) == 3:
        aa = translate_codon(codon)
        if aa == '*':
            break
        protein += aa

print(f"Nucleotide sequence length: {len(nt_sequence)}")
print(f"Protein sequence length: {len(protein)}")
print(f"Protein sequence: {protein[:50]}...")

# Test mutational scan with NT-level mutations
print("\n" + "="*80)
print("Testing nucleotide-level mutation scanning...")
print("="*80)

result = mutational_scan(
    sequence=protein,
    region_start=1,
    region_end=10,  # Scan first 10 amino acids
    gene_hint="rpoB",
    enable_drug_filtering=True,
    resistance_cutoff=0.7,
    demo_mode=True,
    nucleotide_sequence=nt_sequence,
    orf_start=0
)

print(f"\nTotal mutations found: {len(result['matrix'])}")

# Show first 10 mutations
print("\nFirst 10 mutations:")
print("-" * 120)
print(f"{'NT Pos':<8} {'WT NT':<6} {'MUT NT':<7} {'WT Codon':<10} {'MUT Codon':<11} {'AA Pos':<7} {'WT AA':<6} {'MUT AA':<7} {'Group':<8} {'Max P(resist)':<15}")
print("-" * 120)

for i, mut in enumerate(result['matrix'][:10]):
    max_p = max(mut['p_resist'].values())
    print(f"{mut['nt_pos']:<8} {mut['wt_nt']:<6} {mut['mut_nt']:<7} {mut['wt_codon']:<10} {mut['mut_codon']:<11} {mut['pos']:<7} {mut['wt']:<6} {mut['mut']:<7} {mut['group']:<8} {max_p:<15.4f}")

print("\n" + "="*80)
print("Test completed successfully!")
print("="*80)

# Test backward compatibility (no NT sequence)
print("\nTesting backward compatibility (AA-level mutations)...")
result_aa = mutational_scan(
    sequence=protein,
    region_start=1,
    region_end=5,
    gene_hint="rpoB",
    enable_drug_filtering=True,
    resistance_cutoff=0.7,
    demo_mode=True
)

print(f"Total mutations (AA-level): {len(result_aa['matrix'])}")
print("Sample mutation:", result_aa['matrix'][0])
