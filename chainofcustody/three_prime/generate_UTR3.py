def generate_mrna_sponge_utr(mirna_sequences, num_sites=16):
    """
    Generates a 3'UTR mRNA sequence with alternating, bulged miRNA sponge sites.
    
    Parameters:
    mirna_sequences (list of str): List of mature miRNA sequences (5' to 3') in A, U, G, C.
                                   (Can also accept a single string for backward compatibility).
    num_sites (int): Total number of sponge sites desired. Default is 16.
    
    Returns:
    dict: Contains the final 3'UTR sequence and the blueprints for each single site.
    """
    # Allow a single string to be passed by converting it to a list
    if isinstance(mirna_sequences, str):
        mirna_sequences = [mirna_sequences]
        
    # Helper 1: Reverse Complement for RNA
    def reverse_complement(seq):
        complement = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
        return ''.join(complement[base] for base in reversed(seq))
        
    # Helper 2: Create a guaranteed mismatch for the bulge
    def create_mismatch(seq):
        mismatch_map = {'A': 'C', 'U': 'G', 'G': 'U', 'C': 'A'}
        return ''.join(mismatch_map[base] for base in seq)
        
    # Process all input miRNAs and store their bulged site blueprints
    sponge_sites = []
    for mirna_seq in mirna_sequences:
        mirna = mirna_seq.upper().replace('T', 'U')
        rc_mirna = reverse_complement(mirna)
        
        # Slice into domains
        seed_match = rc_mirna[-8:]
        bulge_rc = rc_mirna[-12:-8]
        three_prime_match = rc_mirna[:-12]
        
        # Mutate the bulge sequence
        bulge_mismatch = create_mismatch(bulge_rc)
        
        # Assemble and store the single bulged site (kept uppercase)
        site = three_prime_match + bulge_mismatch + seed_match
        sponge_sites.append(site)
    
    # Define non-homologous, low-structure spacers (all lowercase)
    spacers = [
        'aauu', 'ucga', 'caag', 'auac', 'gaau',
        'cuua', 'uuca', 'agcu', 'uacg', 'gaua',
        'cuac', 'acuc', 'uguu', 'caua', 'ucuu', 'agau'
    ]
    
    # Assemble the multi-site cassette
    cassette = ""
    for i in range(num_sites):
        # Alternate through the generated sponge sites
        current_site = sponge_sites[i % len(sponge_sites)]
        cassette += current_site
        
        # Add a spacer after every site except the last one
        if i < num_sites - 1:
            cassette += spacers[i % len(spacers)]
            
    # Build the final 3'UTR environment with your custom sequence
    stop_codon = "UAA"
    lead_in = "gcauac"   
    lead_out = "gauc"    
    poly_a_signal = "CUCAGGUGCAGGCUGCCUAUCAGAAGGUGGUGGCUGGUGUGGCCAAUGCCCUGGCUCACAAAUACCACUGAGAUCUUUUUCCCUCUGCCAAAAAUUAUGGGGACAUCAUGAAGCCCCUUGAGCAUCUGACUUCUGGCUAAUAAAGGAAAUUUAUUUUCAUUGCAAUAGUGUGUUGGAAUUUUUUGUGUCUCUCACUCGGAAGGACAUAUGGGAGGGCAAAUCAUUUAAAACAUCAGAAUGAGUAUUUGGUUUAGAGUUUGGCA"
    
    final_utr = f"{stop_codon}{lead_in}{cassette}{lead_out}{poly_a_signal}"
    
    return {
        "single_sites": sponge_sites,
        "full_utr": final_utr
    }
