#!/usr/bin/env python3

# Import libraries
import re
import sys

# Stop codons list
stop_codon_lst = ["TAA", "TAG", "TGA"]

# Codon table
codon_dict = {
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
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
    }

# Transfer the user input file to the fname variable. If there is none, terminate the program
if len(sys.argv) == 2:
	fname = sys.argv[1]
elif len(sys.argv) == 1:
	print("ERROR: Not input file was selected. Please, select a valid FASTA file")
	sys.exit()
elif len(sys.argv) > 2:
	print("ERROR? This program only accepts one input.", len(sys.argv) - 1, "were provided.")
	sys.exit()

# Check if input file exists. If not, terminate the program
while True:
	try:
		fhand = open(fname)
		fhand_check = open(fname) # This second opened file will be used for line counting and checking for "non-nucleotide" characters
		break
	except:
		print("ERROR: Input file name does not exist. Try again.")
		sys.exit()

# Check the number of lines in the input file (useful for extracting the last sequence based on the proposed logic). In addition, check if there are any "non" nucleotide characters)
line_count_fixed = 0
non_nt_count = 0
for lines in fhand_check:
	lines = lines.strip()	
	line_count_fixed += 1
	if not lines.startswith(">"):
		non_nt = re.findall("[^ATCG]", lines) # Matched any "non-nucleotide" character
		if len(non_nt) > 0:
			non_nt_count += 1

#  Some more variables and lists
line_count = 0 # It will be used for comparison againsT line_count_fixed in order to extract the last sequence
ID_lst = list() # Append the found IDs to it. It will be required for printing the ID line before the sequence. Furthermore, it will take part in the logic to get new sequences.
seq = None # Variable intended to gather the nucleotide sequences in each of the valid file's lines 
seq_lst = list() # It will be used for appending the sequences gathered by the "seq" variable. Doing this, it's possible obtain the reverse complements. Furthermore, it will be required to get the sequences' frames.
rev_seq_lst = list() # Append the reverse complements into this list

# Parse the FASTA file a second time to extract the longest ORFs
for lines in fhand:
	if non_nt_count > 0: # It avoids potential errors in the presence of "non-nucleotide" characters by terminating the program.
		print("Non nucleotide character was found in at least one of your sequences. Terminating process.")
		sys.exit()
	
	# Iterate through the lines to get the IDs and sequences
	lines = lines.strip()
	line_count += 1
	lines = lines.upper()
	if lines.startswith(">"):
		ID_lst.append(lines)
	else: # Collect the sequences
		if seq is None:
			seq = lines
		else:
			seq += lines

	if len(ID_lst) > 1 or line_count == line_count_fixed: # Logic used to get individual genes and the last sequences without errors
		seq_lst.append(seq)

		# Generate the reverse complement sequence and append it to a list
		seq_reverse_comp = seq
		seq_reverse_comp = seq_reverse_comp.replace("A", "t")
		seq_reverse_comp = seq_reverse_comp.replace("T", "a")
		seq_reverse_comp = seq_reverse_comp.replace("C", "g")
		seq_reverse_comp = seq_reverse_comp.replace("G", "c")
		seq_reverse_comp = seq_reverse_comp[::-1]
		seq_reverse_comp = seq_reverse_comp.upper()
		rev_seq_lst.append(seq_reverse_comp)
		
		# Determine the 6 possible frames
		# (+) strand frames
		frame_1 = seq_lst[0][0:]
		frame_2 = seq_lst[0][1:]
		frame_3 = seq_lst[0][2:]

		# (-) strand frames
		frame_4 = rev_seq_lst[0][0:]
		frame_5 = rev_seq_lst[0][1:]
		frame_6 = rev_seq_lst[0][2:]

		# Create a list of frames to iterate through it later
		frames = [frame_1, frame_2, frame_3, frame_4, frame_5, frame_6]
		seq_lst.clear() # Clear the list for new sequences
		rev_seq_lst.clear() # Clear the list for new sequences
		frame_count = 0 # Reset frame_count for other genes to come
		
		# Parse each one of the 6 possible frames and look for ORFs
		ORF_lst = list()
		ORF_dict = {}
		for frame in frames:
			start_pos = 0
			frame_count += 1				
			codons = re.findall("...", frame) # Separate the obtained sequence into codons

			# Iterate through the codons list and get the sequences starting with ATG
			for codon in codons:
				start_pos += 1
				if codon == "ATG":
					ATGs_lst = codons[start_pos-1:]

					# Iterate through the elements of the newly found sequence starting with ATG
					for element in ATGs_lst:									
						if element in stop_codon_lst: # Indicates the obtention of an ORF
							ORF_lst.append(element) # Append the stop codon to the ORF

							# Create a dictionary to store the ORFs and their data (i.e., frame number and coordinates) 
							# I found some contrasting information about the best way to represent the coordinates. I ended using the same as in https://www.bioinformatics.org/sms2/orf_find.html, since it appeared more logic to me
							if "".join(ORF_lst) not in ORF_dict: # Each one of the found ORFs will be unique. Thus, they make great unique keys for a dictionary
								if frame_count == 1:
									ORF_dict["".join(ORF_lst)] = {"FRAME": frame_count, "START": (start_pos * 3) - 2, "END": (len(ORF_lst) * 3) + (start_pos * 3) - 3}
								elif frame_count == 2:
									ORF_dict["".join(ORF_lst)] = {"FRAME": frame_count, "START": (start_pos * 3) - 1, "END": (len(ORF_lst) * 3) + (start_pos * 3) - 2}
								elif frame_count == 3:
									ORF_dict["".join(ORF_lst)] = {"FRAME": frame_count, "START": (start_pos * 3), "END": (len(ORF_lst) * 3) + (start_pos * 3) - 1}
								elif frame_count == 4:
									ORF_dict["".join(ORF_lst)] = {"FRAME": frame_count, "START": (start_pos * 3) - 2, "END": (len(ORF_lst) * 3) + (start_pos * 3) - 3}
								elif frame_count == 5:
									ORF_dict["".join(ORF_lst)] = {"FRAME": frame_count, "START": (start_pos * 3) - 1, "END": (len(ORF_lst) * 3) + (start_pos * 3) - 2}
								elif frame_count == 6:
									ORF_dict["".join(ORF_lst)] = {"FRAME": frame_count, "START": (start_pos * 3), "END": (len(ORF_lst) * 3) + (start_pos * 3) - 1}									
								
								# Clear the below lists for new sequences
								ORF_lst.clear()
								ATGs_lst.clear()
						else:
							ORF_lst.append(element) # If no stop codon is found, continue appending codons to the putative ORF

		# Extract the longest ORF and its attributes (i.e., frame number and coordinates)
		longest_ORF = None
		for key, value in ORF_dict.items():
			if longest_ORF is None or len(key) > len(longest_ORF):
				longest_ORF = key
				frame = value["FRAME"]
				start_pos = value["START"]
				end_pos = value["END"]
		
		# Save the longest ORF nucleotide sequence
		with open("ORF.fna", "a") as f:
			print(ID_lst[0] + "_frame" + str(frame) + "_" + str(start_pos) + "_" + str(end_pos), file=f)
			print(longest_ORF, file=f) 
					
		# Translate the longest ORF's nucleotide sequence 
		aa_sequence = None
		codons = re.findall("...", longest_ORF) # Separate the ORF into codons
		for codon in codons:
			if codon in codon_dict: # Check the codon correspondence in the codon_dict
				aa = codon_dict[codon]
				if aa_sequence is None:
					aa_sequence = aa
				else:
					aa_sequence = aa_sequence + aa

		# Save the longest ORF's translated sequence
		with open("ORF.faa", "a") as f:
			print(ID_lst[0] + "_frame" + str(frame) + "_" + str(start_pos) + "_" + str(end_pos), file=f)
			print(aa_sequence, file=f)
		
		# Clear all the remaining variables and lists for new sequences
		ID_lst.pop(0)
		ORF_dict.clear()
		seq = None
		longest_ORF = None
		aa_sequence = None
