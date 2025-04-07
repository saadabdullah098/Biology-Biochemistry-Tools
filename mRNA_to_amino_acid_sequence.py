"""
    This program takes a user inpute of mRNA and converts it to a amino acid sequence.
    To run the program, run this script on the terminal or any code editor. 
    
    Note: The program assumes that the first AUG is the start codon and begins translating from there and stops at the first stop codon. 
    To further modify the program, add in a secondary verification of checking for the Kozak sequence of the organism the protein originates from to ensure the correct location is selected
    for the start codon. 
"""



def is_valid_mrna(mrna: str) -> bool:
    """
    Check if a string is a valid mRNA sequence (contains only A, U, G, C).
    Returns True if valid, False otherwise.
    """
    # Convert to uppercase first
    mrna = mrna.upper()
    # Check if string contains only A, U, G, C
    return all(base in 'AUGC' for base in mrna)

def translate_mrna(mrna_sequence):
    """
    Translates an mRNA sequence into a chain of amino acids.
    
    Args:
        mrna_sequence (str): The mRNA sequence to translate
        
    Returns:
        list: A list of amino acids in the translated protein
    """
    # Define the codon table
    CODON_TABLE = {
        'AUG': 'Met', 'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
        'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
        'UAU': 'Tyr', 'UAC': 'Tyr', 'UGU': 'Cys', 'UGC': 'Cys', 'UGG': 'Trp',
        'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
        'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
        'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
        'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
        'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile',
        'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
        'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
        'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
        'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
        'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
        'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
        'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly',
        'UAA': 'STOP', 'UAG': 'STOP', 'UGA': 'STOP'
    }
    
    # Find the start codon (first occurrence of AUG)
    start_index = mrna_sequence.find('AUG')
    
    # If no start codon is found, return an empty list
    if start_index == -1:
        return []
    
    # Initialize the list to store amino acids
    protein = []
    
    # Start translating from the first AUG
    i = start_index
    
    # Process the sequence in triplets until the end or a stop codon
    while i <= len(mrna_sequence) - 3:
        # Get the current codon
        codon = mrna_sequence[i:i+3]
        
        # Translate the codon to an amino acid
        amino_acid = CODON_TABLE.get(codon)
        
        # If it's a stop codon, end translation
        if amino_acid == 'STOP':
            break
        
        # Add the amino acid to the protein
        protein.append(amino_acid)
        
        # Move to the next codon
        i += 3
    
    return protein

def main():
    """
    Main function to handle user input for mRNA translation and display the results.
    """
    print("mRNA Sequence Translation Tool")
    print("=============================")
    print("This program translates mRNA sequences into amino acid chains.")
    
    while True:
        print("\nOptions:")
        print("1. Translate an mRNA sequence")
        print("2. Exit program")
        
        choice = input("\nEnter your choice (1-2): ")
        
        if choice == '2':
            print("Exiting program. Goodbye!")
            break
        
        elif choice == '1':
            mrna_input = input("\nEnter an mRNA sequence: ").strip().upper()
            
            if not mrna_input:
                print("Error: Empty input. Please enter a valid mRNA sequence.")
                continue
            
            if not is_valid_mrna(mrna_input):
                print("Error: Invalid mRNA sequence. Sequence must contain only A, U, G, C.")
                continue
            
            protein = translate_mrna(mrna_input)
            
            print("\nTranslation Results:")
            print("-------------------")
            
            if not protein:
                print("No start codon (AUG) found. Could not translate sequence.")
            else:
                print(f"Start codon position: {mrna_input.find('AUG')}")
                print(f"Protein sequence: {' - '.join(protein)}")
                print(f"Number of amino acids: {len(protein)}")
                
                # Find where translation stopped
                start_index = mrna_input.find('AUG')
                end_index = start_index + (len(protein) * 3)
                
                if end_index + 3 <= len(mrna_input):
                    stop_codon = mrna_input[end_index:end_index+3]
                    print(f"Translation stopped at: {stop_codon} (Stop codon)")
                else:
                    print("Translation stopped at: End of sequence (No stop codon found)")
        
        else:
            print("Invalid choice. Please enter 1 or 2.")

if __name__ == "__main__":
    main()
