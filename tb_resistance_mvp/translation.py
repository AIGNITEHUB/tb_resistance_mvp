"""
DNA/RNA to Protein Translation
"""

# Standard genetic code
CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}


def translate_dna(sequence: str, frame: int = 0) -> str:
    """
    Translate DNA sequence to protein
    
    Args:
        sequence: DNA sequence (ATCG)
        frame: Reading frame (0, 1, or 2)
    
    Returns:
        Protein sequence (1-letter amino acids)
    """
    # Clean sequence
    sequence = sequence.upper().replace(" ", "").replace("\n", "")
    sequence = sequence.replace("U", "T")  # Handle RNA input
    
    # Validate
    valid_bases = set("ATCG")
    if not all(base in valid_bases for base in sequence):
        raise ValueError("Invalid DNA sequence. Only A, T, C, G allowed.")
    
    # Apply reading frame
    sequence = sequence[frame:]
    
    # Translate
    protein = []
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        aa = CODON_TABLE.get(codon, 'X')  # X for unknown
        if aa == '*':  # Stop codon
            break
        protein.append(aa)
    
    return ''.join(protein)


def detect_orf(sequence: str, min_length: int = 100) -> list:
    """
    Detect Open Reading Frames (ORFs)
    
    Args:
        sequence: DNA sequence
        min_length: Minimum ORF length in nucleotides
    
    Returns:
        List of ORFs with start, end, frame, protein
    """
    sequence = sequence.upper().replace("U", "T")
    orfs = []
    
    # Check all 3 frames
    for frame in range(3):
        seq_frame = sequence[frame:]
        
        for i in range(0, len(seq_frame) - 2, 3):
            codon = seq_frame[i:i+3]
            
            # Start codon
            if codon == 'ATG':
                start = i + frame
                
                # Find stop codon
                for j in range(i, len(seq_frame) - 2, 3):
                    stop_codon = seq_frame[j:j+3]
                    if stop_codon in ('TAA', 'TAG', 'TGA'):
                        end = j + frame + 3
                        
                        if (end - start) >= min_length:
                            orf_seq = sequence[start:end]
                            protein = translate_dna(orf_seq, frame=0)
                            
                            orfs.append({
                                'start': start + 1,  # 1-indexed
                                'end': end,
                                'frame': frame,
                                'length_nt': end - start,
                                'length_aa': len(protein),
                                'protein': protein
                            })
                        break
    
    return orfs


# Test
if __name__ == "__main__":
    # Test translation
    dna = "ATGGCGTCAACAAAGCAGCTGCTGGCGGTCCATGTGCCGCGTAAC"
    protein = translate_dna(dna)
    print(f"DNA: {dna}")
    print(f"Protein: {protein}")
    print(f"Expected: MASTKQLLAVGHVPRN")
    
    # Test ORF detection
    test_seq = "AAAAAATGGCGTCAACAAAGCAGCTGCTGGCGGTCCATGTGCCGCGTAACTAA"
    orfs = detect_orf(test_seq, min_length=30)
    print(f"\nFound {len(orfs)} ORFs:")
    for orf in orfs:
        print(f"  Frame {orf['frame']}: {orf['start']}-{orf['end']} ({orf['length_aa']} AA)")