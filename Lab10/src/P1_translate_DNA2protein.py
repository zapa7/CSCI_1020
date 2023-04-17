import sys
import optparse

# refer to https://www.geeksforgeeks.org/dna-protein-python-3/
print('')
print('  ############################################################################')
print('  #                                                                          #')
print('  #             Python program to convert DNA to protein                     #')
print('  #                     CSCI 1020                                            #')
print('  #                                                                          #')
print('  #       Example: python P1_translate_DNA2protein.py  --dna DNA1.fasta      #')
print('  #                                                                          #')
print('  ############################################################################')
print('')


#############  (1)  Define and parsing command-line options  ###########  

parser = optparse.OptionParser()
parser.add_option('--dna', dest='dna',
        default = '',   
        help = 'FASTA file containing the DNA sequence')


(options,args) = parser.parse_args()

if not options.dna:
        print('Error ! DNA sequence not defined. Exiting application...')
        sys.exit(1)


#############  (2)  Get path of DNA file ###########  
DNA_file = options.dna

#############  (3)  Define Transcription and translate functions for DNA -> RNA -> Protein ###########  
def transcribe(seq): 
	RNA_seq = seq.replace('T', 'U')
	return RNA_seq 

def translate(seq): 
	
	table = { 
		'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M', 
		'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T', 
		'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K', 
		'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',				 
		'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L', 
		'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P', 
		'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q', 
		'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R', 
		'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V', 
		'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A', 
		'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E', 
		'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G', 
		'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S', 
		'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L', 
		'UAC':'Y', 'UAU':'Y', 'UAA':'_', 'UAG':'_', 
		'UGC':'C', 'UGU':'C', 'UGA':'_', 'UGG':'W', 
	} 
	protein ="" 
	### Define location of start codon (AUG)
	start_codon = seq.find('AUG')
	print("Start codon is found in RNA at %i-th position" % (start_codon+1))
	for i in range(start_codon, len(seq), 3): 
		codon = seq[i:i + 3]
		if table[codon] == '_': ## this is stop codon
			print("Stop codon is found in RNA at %i-th position and ends at %i-th position" % ((i+1),(i+3)))
			break
		protein+= table[codon] 
	return protein 


def read_seq(inputfile): 
	f = open(inputfile, 'r')
	seq = f.readlines()
	f.close()
	# removing the trailing "\n" and any header lines
	seq = [line.strip() for line in seq if not '>' in line]
	seq = ''.join( seq )
	return seq 



#############  (4)  Loading DNA sequence ###########  

DNA_sequence = read_seq(DNA_file) 

print("(I) Loading DNA sequence\n");
print("DNA sequence: %s" % DNA_sequence);
print("\n")
print("Length of DNA sequence is: %i\n\n" % (len(DNA_sequence)))


#############  (5)  Get RNA sequence from DNA ###########

print("(II) Transcribe DNA to RNA sequence\n")
RNA_sequence = transcribe(DNA_sequence) 
print("RNA sequence: %s" % RNA_sequence)
print("\n")
print("Length of RNA sequence is: %i\n\n" % (len(RNA_sequence)))

#############  (6)  Get Protein sequence from RNA ###########  
print("(III) Translate RNA to Protein sequence\n");

Protein_sequence = translate(RNA_sequence) 
print("\n")
print("Protein sequence: %s" % Protein_sequence)
print("\n")
print("Length of Protein sequence is: %i\n\n" % (len(Protein_sequence)))
