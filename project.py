'''
This project is a primer generator for a DNA sequence.
Program will design primers with lengths of 18-24 bases.
Melting temperature (Tm) will then be calculated for each of the primers.
Formula: 4°C*(# G/C nucleotides) + 2°C*(# A/T nucleotides).
Ideally 50 - 60°C.
'''

# Sequence can be pasted into dna.txt file for better access than this code.
with open('dna.txt') as file:
    for line in file:
        RAW_SEQUENCE = file.read()

# This is the DNA sequence for Green Fluorescent Protein (GFP)
# RAW_SEQUENCE = "atgGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAA GTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGC TGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAG CAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTA CAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGG ACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAAC GGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACAC CCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACG AGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAG"


def main():
    dna_sequence = clean_sequence(RAW_SEQUENCE)
    fwd_prime, fwd_tm = generate_fwd_primer(dna_sequence)
    rvr_prime, rvr_tm = generate_rvr_primer(dna_sequence)
    print("For DNA sequence:")
    print(dna_sequence)
    print("Your forward primer is " + str(fwd_prime) + " with length " + str(len(fwd_prime)) + " base pairs,")
#    fwd_tm = melting_temp(fwd_prime)
    print("with melting temperature of " + str(fwd_tm) + "°C")
    print("Your reverse primer is " + str(rvr_prime) + " with length " + str(len(rvr_prime)) + " base pairs,")
#    rvr_tm = melting_temp(rvr_prime)
    print("with melting temperature of " + str(rvr_tm) + "°C")
    print("Your expected PCR product size is " + str(len(dna_sequence)) + " base pairs.")

def clean_sequence(sequence):
    dna_sequence = ""
    for chr in sequence:
        if chr.isalpha():
            dna_sequence = dna_sequence + chr
    dna_sequence = dna_sequence.upper()
    return dna_sequence

def generate_fwd_primer(sequence):
    fwd_primer = ""
    fwd_tm = 0
    for i in range(len(sequence)):
        fwd_primer = fwd_primer + sequence[i]
        if sequence[i] == "G" or sequence[i] == "C":
            fwd_tm += 4
        else:
            fwd_tm += 2
        if fwd_tm > 53:
            break
    return fwd_primer, fwd_tm

def generate_rvr_primer(sequence):
    rvr_sequence = ""
    for i in range(len(sequence)):
        rvr_sequence = sequence[i] + rvr_sequence
    rvr_template = rvr_sequence.replace('A', 'X')
    rvr_template = rvr_template.replace('T', 'A')
    rvr_template = rvr_template.replace('X', 'T')
    rvr_template = rvr_template.replace('G', 'X')
    rvr_template = rvr_template.replace('C', 'G')
    rvr_template = rvr_template.replace('X', 'C')
    rvr_primer = ""
    rvr_tm = 0
    for i in range(len(rvr_template)):
        rvr_primer = rvr_primer + rvr_template[i]
        if rvr_template[i] == "G" or rvr_template[i] == "C":
            rvr_tm += 4
        else:
            rvr_tm += 2
        if rvr_tm > 53:
            break
    return rvr_primer, rvr_tm

'''

No longer needed. Length of primer is determined by melting temp, Tm. Addition of base pairs to primer stops
when Tm is greater than 54°C and melting point is calculated and returned in the generate_primer function.

def melting_temp(primer):
    no_of_GCs = primer.replace('A', '')
    no_of_GCs = no_of_GCs.replace('T', '')
    temp_of_GCs = len(no_of_GCs) * 4
    no_of_ATs = ''
    no_of_ATs = primer.replace('G', '')
    no_of_ATs = no_of_ATs.replace('C', '')
    temp_of_ATs = len(no_of_ATs) * 2
    Tm = temp_of_GCs + temp_of_GCs
    return Tm
'''

if __name__ == '__main__':
    main()