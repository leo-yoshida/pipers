#! /usr/bin/python3
#Usage: chomp.py file.fasta number_slices

import math
import sys

n_fatias = int(sys.argv[2])

def count_fasta(arquivo_fasta):
	flecha = 0
	linhas = 0
	with open(arquivo_fasta, "r") as f:
		for line in f:
			if line.startswith('>'):
				flecha += 1
				linhas += 1
			else:
				linhas += 1
	return flecha, linhas

n_sequencias, n_linhas_total = count_fasta(sys.argv[1])

step = int( math.ceil(n_sequencias/n_fatias) )

count_sequencias = 0
count_linhas = 0
bloco = 1

with open(sys.argv[1], "r") as file:
	for line in file:
		if line.startswith('>'):
			if count_linhas == 0:
				f = open(sys.argv[1] + "_chomp_" + str(bloco) + ".fasta", "w")
				f.write(line)
				count_linhas += 1
				count_sequencias += 1
			elif count_sequencias == step:
				f.close()
				bloco += 1
				f = open(sys.argv[1] + "_chomp_" + str(bloco) + ".fasta", "w")
				f.write(line)
				count_sequencias = 1
				count_linhas += 1
			elif count_linhas != 0:
				f.write(line)
				count_linhas += 1
				count_sequencias += 1
		else:
			f.write(line)
			count_linhas += 1
			if count_linhas == n_linhas_total:
				f.close()

