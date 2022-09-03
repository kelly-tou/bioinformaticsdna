from collections import Counter
from structures import *

def validate(dna):
	upper_dna = dna.upper()
	for nucleotide in upper_dna:
		if nucleotide not in nucleotides:
			return False
	if not len(upper_dna) == 30:
		return False
	return upper_dna

def nucleotide_frequency(dna):
	frequency_dict = {"A": 0, "C": 0, "G": 0, "T": 0}
	for nucleotide in dna:
		frequency_dict[nucleotide] += 1
	return frequency_dict

def transcription(dna):
	return dna.replace("T", "U")

def rev_complement(dna):
	return ''.join([complement[nucleotide] for nucleotide in dna])[::-1]

def gc_content(dna):
	return round((dna.count("C") + dna.count("G")) / len(dna) * 100)

def gc_content_sub(dna, k=5):
	res = []
	for i in range(0, len(dna) - k +1, k):
		sub = dna[i:i+k]
		res.append(gc_content(sub))
	return res

def translation(dna, init_pos = 0):
	return [dna_codons[dna[pos: pos+3]] for pos in range(init_pos, len(dna)-2, 3)]

def codon(dna, aminoacid):
	tmp_list = []
	for i in range(0, len(dna) - 2, 3):
		if dna_codons[dna[i:i+3]] == aminoacid:
			tmp_list.append(dna[i:i+3])
	freq_dict = dict(Counter(tmp_list))
	total = sum(freq_dict.values())
	for dna in freq_dict:
		freq_dict[dna] = round(freq_dict[dna] / total, 2)
	return freq_dict

def reading_frames(dna):
	frames = []
	frames.append(translation(dna, 0))
	frames.append(translation(dna, 1))
	frames.append(translation(dna, 2))
	frames.append(translation(rev_complement(dna), 0))
	frames.append(translation(rev_complement(dna), 1))
	frames.append(translation(rev_complement(dna), 2))
	return frames

def proteins(aa_dna):
	current = []
	proteins = []
	for aa in aa_dna:
		if aa == "_":
			if current:
				for p in current:
					proteins.append(p)
				current = []
		else:
			if aa == "M":
				current.append("")
			for i in range(len(current)):
				current[i] += aa
	return proteins

def all_proteins(dna, start_read_pos = 0, end_read_pos = 0, ordered = False):
	if end_read_pos > start_read_pos:
		rfs = reading_frames(dna[start_read: end_read])
	else:
		rfs = reading_frames(dna)
	res = []
	for rf in rfs:
		pros = proteins(rf)
		for p in pros:
			res.append(p)
	if ordered:
		return sorted(res, key = len, reverse = True)
	return res
