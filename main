import pygame
from dnafunctions import *
from structures import *
import random
from text import Button
import sys

pygame.init()
scr_w = 1200
scr_h = 700
screen = pygame.display.set_mode((scr_w, scr_h))
clock = pygame.time.Clock()
SCREEN_UPDATE = pygame.USEREVENT
pygame.time.set_timer(SCREEN_UPDATE, 500)
pygame.display.set_caption("Bioinformatics: DNA/RNA Sequence Analysis")

print("Press r for randomly generated DNA, or input your own. DNA must be 30 nucleotides long")
rand_dna = str
k = input()
if k == "r":
	rand_dna = ''.join([random.choice(nucleotides) for nucleotide in range(30)])
else:
	rand_dna = k

if validate(rand_dna):
	dna = validate(rand_dna)
else:
	print("Invalid Input. DNA must be 30 nucleotides long")
	exit()

bar = pygame.image.load("nuc_freq.png").convert_alpha()
bar = pygame.transform.scale(bar, (350, 300))

text_size = 40
sequence = Button(scr_w/2 - 20, 0, text_size, text_size, text = f"DNA Sequence: {dna}", text_color = (237,111,248))
length = Button(scr_w/2 - 20, text_size - 5, text_size, text_size, text = f"Length: {len(dna)}", text_color = (237,111,248))

bar_rect = pygame.Rect(10, 75, 50, 50)

transcription_label = Button(scr_w/2 - 66, 130,0,0, text = "DNA/RNA Transcription:", text_color = (103,222,130))
transcription = Button(scr_w/2 + 300, 130,0,0, text = transcription(dna), text_color = (255,255,255))

dna_label = Button(scr_w/2 - 8, 170,0,0, text = "DNA String:", text_color = (103,222,130))
dna_string = Button(scr_w/2 + 300, 170,0,0, text = f"5'  {dna}  3'", text_color = (255,255,255))

connector = Button(scr_w/2 + 305, 190,0,0, text = ''.join([ "| " for c in range(len(dna))]), text_color = (255,255,255))

complement_label = Button(scr_w/2 - 17, 220,0,0, text = "Complement:", text_color = (103,222,130))
complement = Button(scr_w/2 + 300, 220,0,0, text = f"3'  {rev_complement(dna)[::-1]}  5'", text_color = (255,255,255))

rev_complement_label = Button(scr_w/2 - 60, 270,0,0, text = "Reverse Complement:", text_color = (103,222,130))
rev_complement = Button(scr_w/2 + 300, 270,0,0, text = f"5'  {rev_complement(dna)}  3'", text_color = (255,255,255))

gc_perc = Button(190, 455,0,0, text = f"{gc_content(dna)}%", text_color = (255,255,255))
gc_label = Button(190, 515,0,0, text = "GC Content", text_color = (255,255,255))

gc = pygame.image.load("gc_sub.png").convert_alpha()
gc = pygame.transform.scale(gc, (600, 200))
gc_rect = pygame.Rect(350, 360, 50, 50)

gc_sub1 = Button(450, 500-gc_content_sub(dna)[0]*2,0,0, text = f"{gc_content_sub(dna)[0]}%", text_color = (90,131,247))
gc_sub2 = Button(517, 500-gc_content_sub(dna)[1]*2,0,0, text = f"{gc_content_sub(dna)[1]}%", text_color = (90,131,247))
gc_sub3 = Button(597, 500-gc_content_sub(dna)[2]*2,0,0, text = f"{gc_content_sub(dna)[2]}%", text_color = (90,131,247))
gc_sub4 = Button(680, 500-gc_content_sub(dna)[3]*2,0,0, text = f"{gc_content_sub(dna)[3]}%", text_color = (90,131,247))
gc_sub5 = Button(764, 500-gc_content_sub(dna)[4]*2,0,0, text = f"{gc_content_sub(dna)[4]}%", text_color = (90,131,247))
gc_sub6 = Button(844, 500-gc_content_sub(dna)[5]*2,0,0, text = f"{gc_content_sub(dna)[5]}%", text_color = (90,131,247))
gc_sub_label = Button(650, 557,0,0, text = "GC Content in Subsection k = 5", text_color = (255,255,255))

amino_seq_label = Button(200, 620,0,0, text = "Aminoacids Sequence from DNA:", text_color = (144,252,244))
amino_seq = Button(550, 620,0,0, text = str(translation(dna, 0)), text_color = (255,255,255))
codon_freq_label = Button(140, 660,0,0, text = "Codon Frequency (L):", text_color = (144,252,244))
codon_freq = Button(270+len(str(codon(dna, "L"))) *4, 660,0,0, text = str(codon(dna, "L")), text_color = (255,255,255))

read_fram_label = Button(1000, 340,0,0, text = "Reading Frames:", text_color = (144,252,244))
read_fram1 = Button(1050, 370,0,0, text = str(reading_frames(dna)[0]), text_color = (255,255,255))
read_fram2 = Button(1050, 395,0,0, text = str(reading_frames(dna)[1]), text_color = (255,255,255))
read_fram3 = Button(1050, 420,0,0, text = str(reading_frames(dna)[2]), text_color = (255,255,255))
read_fram4 = Button(1050, 445,0,0, text = str(reading_frames(dna)[3]), text_color = (255,255,255))
read_fram5 = Button(1050, 470,0,0, text = str(reading_frames(dna)[4]), text_color = (255,255,255))
read_fram6 = Button(1050, 495,0,0, text = str(reading_frames(dna)[5]), text_color = (255,255,255))

protein_label = Button(1050, 565,0,0, text = "All Proteins in 6 Open", text_color = (144,252,244))
protein_label2 = Button(1050, 595,0,0, text = "Reading Frames:", text_color = (144,252,244))
proteins = []
for pro in all_proteins(dna, 0, 0, True):
	proteins.append(str(pro))
protein = Button(1050, 625,0,0, text = str(proteins), text_color = (255,255,255))

while True:
	screen.fill((22,59,109))
	sequence.draw(screen, 36)
	length.draw(screen, 36)

	screen.blit(bar, bar_rect)
	pygame.draw.rect(screen, (243,175,110), (90,325-nucleotide_frequency(dna)["A"]*15,40,nucleotide_frequency(dna)["A"]*15))
	pygame.draw.rect(screen, (86,188,185), (150,325-nucleotide_frequency(dna)["T"]*15,40,nucleotide_frequency(dna)["T"]*15))
	pygame.draw.rect(screen, (252,238,135), (207,325-nucleotide_frequency(dna)["G"]*15,40,nucleotide_frequency(dna)["G"]*15))
	pygame.draw.rect(screen, (237,112,107), (261,325-nucleotide_frequency(dna)["C"]*15,40,nucleotide_frequency(dna)["C"]*15))

	pygame.draw.rect(screen, (40,78,123), (400,100, 750, 200))
	transcription_label.draw(screen, 30)
	transcription.draw(screen, 30)
	dna_label.draw(screen, 30)
	dna_string.draw(screen, 30)
	connector.draw(screen, 35)
	complement_label.draw(screen, 30)
	complement.draw(screen, 30)
	rev_complement_label.draw(screen, 30)
	rev_complement.draw(screen, 30)

	pygame.draw.rect(screen, (63,98,137), (50,390,275,150))

	gc_perc.draw(screen, 150)
	gc_label.draw(screen, 30)

	screen.blit(gc, gc_rect)
	pygame.draw.rect(screen, (90,131,247), (430,511-gc_content_sub(dna)[0]*2,30,gc_content_sub(dna)[0]*2))
	pygame.draw.rect(screen, (90,131,247), (500,511-gc_content_sub(dna)[1]*2,30,gc_content_sub(dna)[1]*2))
	pygame.draw.rect(screen, (90,131,247), (580,511-gc_content_sub(dna)[2]*2,30,gc_content_sub(dna)[2]*2))
	pygame.draw.rect(screen, (90,131,247), (663,511-gc_content_sub(dna)[3]*2,30,gc_content_sub(dna)[3]*2))
	pygame.draw.rect(screen, (90,131,247), (745,511-gc_content_sub(dna)[4]*2,30,gc_content_sub(dna)[4]*2))
	pygame.draw.rect(screen, (90,131,247), (827,511-gc_content_sub(dna)[5]*2,30,gc_content_sub(dna)[5]*2))

	gc_sub1.draw(screen, 30)
	gc_sub2.draw(screen, 30)
	gc_sub3.draw(screen, 30)
	gc_sub4.draw(screen, 30)
	gc_sub5.draw(screen, 30)
	gc_sub6.draw(screen, 30)
	gc_sub_label.draw(screen, 30)

	amino_seq_label.draw(screen, 30)
	codon_freq_label.draw(screen, 30)
	amino_seq.draw(screen, 30)
	codon_freq.draw(screen, 30)
	read_fram_label.draw(screen, 30)
	read_fram1.draw(screen, 22)
	read_fram2.draw(screen, 22)
	read_fram3.draw(screen, 22)
	read_fram4.draw(screen, 22)
	read_fram5.draw(screen, 22)
	read_fram6.draw(screen, 22)

	pygame.draw.rect(screen, (86,116,152), (923,545, 250,100))
	protein_label.draw(screen, 30)
	protein_label2.draw(screen, 30)
	protein.draw(screen, 30)

	for e in pygame.event.get():
		if e.type == pygame.QUIT:
			pygame.quit()
			sys.exit()

	pygame.display.update()
