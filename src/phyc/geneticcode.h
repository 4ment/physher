// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef GENETIC_CODE_H
#define GENETIC_CODE_H

// Taken from BEAST GeneticCode.java

static const char * const GENETIC_CODE_TABLES[15] = {
	// Universal
	"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF",
	// Vertebrate Mitochondrial
	"KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
	// Yeast
	"KNKNTTTTRSRSMIMIQHQHPPPPRRRRTTTTEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
	// Mold Protozoan Mitochondrial
	"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
	// Mycoplasma
	"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
	// Invertebrate Mitochondrial
	"KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
	// Ciliate
	"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSS*CWCLFLF",
	// Echinoderm Mitochondrial
	"NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
	// Euplotid Nuclear
	"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSCCWCLFLF",
	// Bacterial
	"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF",
	// Alternative Yeast
	"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLSLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF",
	// Ascidian Mitochondrial
	"KNKNTTTTGSGSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
	// Flatworm Mitochondrial
	"NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYY*YSSSSWCWCLFLF",
	// Blepharisma Nuclear
	"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YQYSSSS*CWCLFLF",
	// No stops
	"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYYQYSSSSWCWCLFLF"
};

static const char * const GENETIC_CODE_NAMES[15] = {
	"Universal",
	"Vertebrate Mitochondrial",
	"Yeast",
	"Mold Protozoan Mitochondrial",
	"Mycoplasma",
	"Invertebrate Mitochondrial",
	"Ciliate",
	"Echinoderm Mitochondrial",
	"Euplotid Nuclear",
	"Bacterial",
	"Alternative Yeast",
	"Ascidian Mitochondrial",
	"Flatworm Mitochondrial",
	"Blepharisma Nuclear",
	"No stops"
};

static int const NUMBER_OF_CODONS[15] = {
	61,
	60,
	62,
	62,
	62,
	62,
	63,
	62,
	62,
	61,
	61,
	62,
	63,
	62,
	64
};

static  char const CODON_TRIPLETS[66][4] = {
	"AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT",
	"AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT",
	"CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT",
	"CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT",
	"GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT",
	"GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT",
	"TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT",
	"TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT",
	"???", "---"
};

void GenticCode_encoding_to_codon_string( int encoding, int genetic_code, char *codon );

#endif
