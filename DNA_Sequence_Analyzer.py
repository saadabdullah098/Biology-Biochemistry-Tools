"""
    This program takes a DNA input and analyzes the base counts and GC/AT contents. 
    To use it, install the dependacies using pip install pandas, collections, re on the terminal. Then run this script and input DNA sequence. 
    
    Dependencies: pandas, collections, re 
"""

import pandas as pd
import collections
import re

def is_valid_dna(dna: str) -> bool:
    """
    Check if a string is a valid DNA sequence (contains only A, T, G, C).
    Returns True if valid, False otherwise.
    """
    # Convert to uppercase first
    dna = dna.upper()
    # Use regex to check if string contains only A, T, G, C
    return bool(re.match(r'^[ATGC]+$', dna))

def dna_makeup(dna: str) -> pd.DataFrame:
    """
    This function takes a DNA sequence input and analyzes the GC content & individual base counts.
    Returns a pandas DataFrame with the analysis results.
    """
    # Validate the DNA sequence
    if not is_valid_dna(dna):
        raise ValueError("Invalid DNA sequence. Sequence must contain only A, T, G, C.")
    
    # Convert to uppercase once
    dna = dna.upper()
    
    # Use Counter for efficient base counting
    base_counts = collections.Counter(dna)
    
    # Get individual base counts
    G = base_counts.get('G', 0)
    C = base_counts.get('C', 0)
    A = base_counts.get('A', 0)
    T = base_counts.get('T', 0)
    
    dna_len = len(dna)
    
    # Calculate percentages
    G_percentage = (G/dna_len) * 100
    C_percentage = (C/dna_len) * 100
    A_percentage = (A/dna_len) * 100
    T_percentage = (T/dna_len) * 100
    GC_content = ((G+C)/dna_len) * 100
    AT_content = ((A+T)/dna_len) * 100
    
    # Create DataFrame with reoriented structure
    data = {
        'Base Count': [A, T, G, C, A+T, G+C],
        'Percentage': [A_percentage, T_percentage, G_percentage, C_percentage, AT_content, GC_content]
    }
    
    # Create index labels
    index_labels = ['A', 'T', 'G', 'C', 'AT Content', 'GC Content']
    
    # Create DataFrame with the reoriented structure
    df = pd.DataFrame(data, index=index_labels)
    
    # Round percentage values to 2 decimal places for better readability
    df['Percentage'] = df['Percentage'].round(2)
    
    return df

def main():
    """
    Main function to handle user input and display the results.
    """
    print("DNA Sequence Analysis Tool")
    print("=========================")
    
    while True:
        dna_input = input("\nEnter a DNA sequence (or 'q' to quit): ")
        
        if dna_input.lower() == 'q':
            print("Exiting program. Goodbye!")
            break
        
        try:
            # Analyze the DNA
            result_df = dna_makeup(dna_input)
            
            # Display the results
            print("\nAnalysis Results:")
            print("----------------")
            print(result_df)
            print("\nTotal sequence length:", len(dna_input))
            
        except ValueError as e:
            print(f"Error: {e}")
            print("Please enter a valid DNA sequence containing only A, T, G, C.")

if __name__ == "__main__":
    main()
