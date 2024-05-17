 # A package designed for DNA, RNA and protein sequence analysis.
 ### Classes provided are Sequence and Invalid. Please import these as - from BioProt import Sequence, Invalid
 ## Fuctions included in the class:
 #### * init method : Raises errors if a sequence is not DNA/RNA sequence but was provided to be so *
   #### - inputs = sequence: String, is_rna: boolean -
 #### * valid_DNAseq() : Checks if the string input is a DNA sequence *
   #### - inputs = NIL -
 #### ~ valid_RNAseq() : Checks if the string input is an RNA sequence
   #### - inputs = NIL -
 #### ~ revseq() : Reverses any sequence
   #### - inputs = NIL -
 #### ~ complseq() : Complements any sequence (DNA/RNA) based on is_rna value
   #### - inputs = NIL -
#### ~ complrevseq() : Reverses and complements any sequence (DNA/RNA)
   #### - inputs = NIL -
#### ~ seqlen() : Returns length of any sequence
   #### - inputs = NIL -
#### ~ countnucs() : Counts the number of nucleotides in any sequence
   #### - inputs = NIL -
#### ~ gc_content() : Gives the GC content as % of the entire sequence
   #### - inputs = NIL -
#### ~ transDNAtoRNA() : Transcribes DNA sequence to RNA
   #### - inputs = NIL -
#### ~ transDNAtoP() : Transforms DNA to Proteins through transcription and translation. Returns only first reading frame (starting from 0)
   #### - inputs = NIL -
#### ~ transDNAtoP_all() : Transforms DNA to Proteins through transcription and translation. Returns all 3 reading frames 
   #### - inputs = NIL -
#### ~ transRNAtoP() : Translates RNA to Protein sequence. Intermediate function for transDNAtoP() and transDNAtoP_all()
   #### - inputs = NIL -
#### ~ pMap(letter) : Returns the protein name corresponding to the letter for a particular codon
  #### - inputs = letter: char -
    
