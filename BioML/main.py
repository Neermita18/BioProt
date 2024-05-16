import string
class Invalid(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)
    
class Sequence:
    dnanucs= {'A', 'C', 'G', 'T'}
    rnanucs= {'A', 'C', 'G', 'U'}
    protein_map = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'AGU': 'S', 'AGC': 'S', 'CCU': 'P', 'CCC': 'P',
        'CCA': 'P', 'CCG': 'P', 'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCU': 'A', 'GCC': 'A',
        'GCA': 'A', 'GCG': 'A', 'UAU': 'Y', 'UAC': 'Y', 'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'UGU': 'C', 'UGC': 'C', 'UGG': 'W', 'CGU': 'R',
        'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G',
        'GGG': 'G'
    }
    def __init__(self, sequence, is_rna=False):
        self.is_rna = is_rna
        if is_rna:
            if not self.valid_RNAseq(sequence):
                raise Invalid("Invalid RNA sequence!")
        else:
            if not self.valid_DNAseq(sequence):
                raise Invalid("Invalid DNA sequence!")
        self.sequence = sequence
        
    def valid_DNAseq(self, sequence):
        for x in sequence:
            if x not in Sequence.dnanucs:
                return False
            
        return True
    def valid_RNAseq(self, sequence):
        for x in sequence:
            if x not in Sequence.rnanucs:
                return False
            
        return True
    def revseq(self):
        return self.sequence[::-1]
    def complseq(self):
        if self.is_rna:
            complement_map = str.maketrans('ACGU', 'UGCA')
        else:
            complement_map = str.maketrans('ACGT', 'TGCA')
        return self.sequence.translate(complement_map)
    def complrevseq(self):
        return self.complseq()[::-1]
            
    
    def seqlen(self):
        return len(self.sequence)

    def countnucs(self):
        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'U':0}
        for nucleotide in self.sequence:
            if nucleotide in counts:
                counts[nucleotide] += 1
            else:
                raise Invalid(f"Invalid nucleotide '{nucleotide}' found in sequence.")
        return counts
    
    def gc_content(self):
        counts = self.countnucs()
        gc_count = counts['G'] + counts['C']
        total_count = len(self.sequence)
        return (gc_count / total_count) * 100 if total_count > 0 else 0

    def transDNAtoRNA(self):
        if self.is_rna:
            raise ValueError("Sequence is already an RNA sequence.")
        return Sequence(self.sequence.replace('T', 'U'), is_rna=True)

    def translate_dna_to_protein(self):
        if self.is_rna:
            raise Invalid("Translation to protein can only be done from a DNA sequence.")
        rna_sequence = self.sequence.replace('T', 'U')
        rna_seq_instance = Sequence(rna_sequence, is_rna=True)
        return rna_seq_instance.transRNAtoP()

    def translate_dna_all_frames(self):
        if self.is_rna:
            raise Invalid("Translation to protein can only be done from a DNA sequence.")
        rna_sequence = self.sequence.replace('T', 'U')
        return [Sequence(rna_sequence[i:], is_rna=True).transRNAtoP() for i in range(3)]

    def transRNAtoP(self):
        protein = []
        for i in range(0, len(self.sequence), 3):
            codon = self.sequence[i:i+3]
            if len(codon) < 3:
                break
            amino_acid = Sequence.protein_map.get(codon)
            if amino_acid == 'Stop':
                break
            if amino_acid:
                protein.append(amino_acid)
        return ''.join(protein)

    def __str__(self):
        return self.sequence
    
dna_seq = Sequence("ATCCCCCCCCCCAAAATGATAGGATGC")
print("Length:", dna_seq.seqlen())
print("Nucleotide Counts:", dna_seq.countnucs())
print("GC Content:", dna_seq.gc_content())
rna_seq = dna_seq.transDNAtoRNA()
print("Transcribed RNA:", rna_seq)
rna_seq = Sequence("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA", is_rna= True)
print("Protein sequence:", rna_seq.transRNAtoP())
dna_seq = Sequence("ATGGCCATGGCGCCCAGAACCTGAGATCAATAGTACCCGTATTAACGGGTGA")
    
print("DNA to protein (start from 0):", dna_seq.translate_dna_to_protein())
print("DNA all reading frames:", dna_seq.translate_dna_all_frames())