import matplotlib.pyplot as plt
import matplotlib.patches as patches
import re

# ---  The 3'UTR Generator Function ---
def generate_mrna_sponge_utr(mirna_sequences, num_sites=16):
    """Generates a 3'UTR mRNA sequence with alternating, bulged miRNA sponge sites."""
    if isinstance(mirna_sequences, str):
        mirna_sequences = [mirna_sequences]
        
    def reverse_complement(seq):
        complement = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
        return ''.join(complement[base] for base in reversed(seq))
        
    def create_mismatch(seq):
        mismatch_map = {'A': 'C', 'U': 'G', 'G': 'U', 'C': 'A'}
        return ''.join(mismatch_map[base] for base in seq)
        
    sponge_sites = []
    for mirna_seq in mirna_sequences:
        mirna = mirna_seq.upper().replace('T', 'U')
        rc_mirna = reverse_complement(mirna)
        
        seed_match = rc_mirna[-8:]
        bulge_rc = rc_mirna[-12:-8]
        three_prime_match = rc_mirna[:-12]
        
        bulge_mismatch = create_mismatch(bulge_rc)
        site = three_prime_match + bulge_mismatch + seed_match
        sponge_sites.append(site)
    
    spacers = [
        'aauu', 'ucga', 'caag', 'auac', 'gaau',
        'cuua', 'uuca', 'agcu', 'uacg', 'gaua',
        'cuac', 'acuc', 'uguu', 'caua', 'ucuu', 'agau'
    ]
    
    cassette = ""
    for i in range(num_sites):
        current_site = sponge_sites[i % len(sponge_sites)]
        cassette += current_site
        if i < num_sites - 1:
            cassette += spacers[i % len(spacers)]
            
    # Note:  synthetic UTRs 
    # often include a backup STOP codon. We will keep the UAA here.
    stop_codon = "UAA"
    lead_in = "gcauac"   
    lead_out = "gauc"    
    poly_a_signal = "CUCAGGUGCAGGCUGCCUAUCAGAAGGUGGUGGCUGGUGUGGCCAAUGCCCUGGCUCACAAAUACCACUGAGAUCUUUUUCCCUCUGCCAAAAAUUAUGGGGACAUCAUGAAGCCCCUUGAGCAUCUGACUUCUGGCUAAUAAAGGAAAUUUAUUUUCAUUGCAAUAGUGUGUUGGAAUUUUUUGUGUCUCUCACUCGGAAGGACAUAUGGGAGGGCAAAUCAUUUAAAACAUCAGAAUGAGUAUUUGGUUUAGAGUUUGGCA"
    
    final_utr = f"{stop_codon}{lead_in}{cassette}{lead_out}{poly_a_signal}"
    return final_utr

# ---  The Plotting Function ---
def plot_mrna_construct(seq_5utr, seq_cds, seq_3utr):
    """Generates a linear map of the mRNA construct directly from sequence strings."""
    
    # Calculate lengths and cumulative positions
    len_5 = len(seq_5utr)
    len_cds = len(seq_cds)
    len_3 = len(seq_3utr)
    total_len = len_5 + len_cds + len_3
    
    start_cds = len_5
    end_cds = start_cds + len_cds
    
    # Setup the Plot
    fig, ax = plt.subplots(figsize=(18, 5))
    ax.set_xlim(-total_len * 0.05, total_len * 1.05)
    ax.set_ylim(0, 1)
    ax.axis('off') 
    
    # Draw the main backbone line
    ax.plot([0, total_len], [0.4, 0.4], color='black', linewidth=2)
    
    # Helper to draw domains
    def draw_domain(start, width, color, label, y=0.3, height=0.2):
        rect = patches.Rectangle((start, y), width, height, facecolor=color, edgecolor='black', zorder=2)
        ax.add_patch(rect)
        ax.text(start + width/2, y + height + 0.05, label, ha='center', va='bottom', fontsize=12, fontweight='bold')

    # Draw Main Domains
    draw_domain(0, len_5, 'lightblue', "5' UTR")
    draw_domain(start_cds, len_cds, 'lightgreen', "CDS")
    draw_domain(end_cds, len_3, 'lightcoral', "3' UTR")
    
    # Annotate Start Codon (First 3 nt of CDS)
    ax.plot([start_cds, start_cds], [0.2, 0.3], color='green', lw=2)
    ax.text(start_cds + 1.5, 0.15, "Start\nCodon", ha='center', va='top', color='green', fontsize=10)
    
    # Annotate Stop Codon (Last 3 nt of CDS)
    ax.plot([end_cds - 3, end_cds - 3], [0.2, 0.3], color='red', lw=2)
    ax.text(end_cds - 1.5, 0.15, "Stop\nCodon", ha='center', va='top', color='red', fontsize=10)

    # Map the miRNA Sponge Sites exclusively in the 3'UTR string
    site_counter = 1
    for match in re.finditer(r'[A-Z]{15,30}', seq_3utr):
        site_start = end_cds + match.start()
        site_width = len(match.group())
        
        # Draw a yellow block over the 3'UTR for the specific site
        rect = patches.Rectangle((site_start, 0.3), site_width, 0.2, facecolor='gold', edgecolor='black', zorder=3)
        ax.add_patch(rect)
        
        # Stagger text heights
        y_text = 0.65 if site_counter % 2 != 0 else 0.85
        ax.text(site_start + site_width/2, y_text, f"Site {site_counter}", ha='center', va='bottom', fontsize=8, rotation=45)
        
        # Draw connecting line
        ax.plot([site_start + site_width/2, site_start + site_width/2], [0.5, y_text - 0.02], color='gray', lw=0.5, zorder=1)
        
        site_counter += 1

    plt.title(f"Synthetic miRNA Sponge Construct ({total_len} nucleotides)", fontsize=16, pad=20)
    plt.tight_layout()
    plt.show()

# ---  Execution Block ---

# Default sequences provided (converted T -> U for RNA consistency)
default_5utr = "GAGTAGTCCCTTCGCAAGCCCTCATTTCACCAGGCCCCCGGCTTGGGGCGCCTTCCTTCCCC".replace('T', 'U')

default_cds = "ATGGCGGGACACCTGGCTTCGGATTTCGCCTTCTCGCCCCCTCCAGGTGGTGGAGGTGATGGGCCAGGGGGGCCGGAGCCGGGCTGGGTTGATCCTCGGACCTGGCTAAGCTTCCAAGGCCCTCCTGGAGGGCCAGGAATCGGGCCGGGGGTTGGGCCAGGCTCTGAGGTGTGGGGGATTCCCCCATGCCCCCCGCCGTATGAGTTCTGTGGGGGGATGGCGTACTGTGGGCCCCAGGTTGGAGTGGGGCTAGTGCCCCAAGGCGGCTTGGAGACCTCTCAGCCTGAGGGCGAAGCAGGAGTCGGGGTGGAGAGCAACTCCGATGGGGCCTCCCCGGAGCCCTGCACCGTCACCCCTGGTGCCGTGAAGCTGGAGAAGGAGAAGCTGGAGCAAAACCCGGAGGAGTCCCAGGACATCAAAGCTCTGCAGAAAGAACTCGAGCAATTTGCCAAGCTCCTGAAGCAGAAGAGGATCACCCTGGGATATACACAGGCCGATGTGGGGCTCACCCTGGGGGTTCTATTTGGGAAGGTATTCAGCCAAACGACCATCTGCCGCTTTGAGGCTCTGCAGCTTAGCTTCAAGAACATGTGTAAGCTGCGGCCCTTGCTGCAGAAGTGGGTGGAGGAAGCTGACAACAATGAAAATCTTCAGGAGATATGCAAAGCAGAAACCCTCGTGCAGGCCCGAAAGAGAAAGCGAACCAGTATCGAGAACCGAGTGAGAGGCAACCTGGAGAATTTGTTCCTGCAGTGCCCGAAACCCACACTGCAGCAGATCAGCCACATCGCCCAGCAGCTTGGGCTCGAGAAGGATGTGGTCCGAGTGTGGTTCTGTAACCGGCGCCAGAAGGGCAAGCGATCAAGCAGCGACTATGCACAACGAGAGGATTTTGAGGCTGCTGGGTCTCCTTTCTCAGGGGGACCAGTGTCCTTTCCTCTGGCCCCAGGGCCCCATTTTGGTACCCCAGGCTATGGGAGCCCTCACTTCACTGCACTGTACTCCTCGGTCCCTTTCCCTGAGGGGGAAGCCTTTCCCCCTGTCTCCGTCACCACTCTGGGCTCTCCCATGCATTCAAACTGA".replace('T', 'U')

# Target miRNA for testing (miR-122-3p)
target_mirna = "AACGCCAUUAUCACACUAAAUA"




# Generate the 3'UTR (with 16 sites)
generated_3utr = generate_mrna_sponge_utr(target_mirna, num_sites=16)

# Generate the plot
plot_mrna_construct(default_5utr, default_cds, generated_3utr)