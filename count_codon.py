#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import division
import sys


def count_codons(seq, codons):
    count = [0 for i in range(64)]
    i = 0
    while (i+1 < len(seq)):
        if seq[i:(i+3)] in codons:
            count[codons[seq[i:(i+3)]]] += 1
        i += 3
    freq = [str(i/(len(seq)/3)) for i in count]
    return freq


seq_file = open(sys.argv[1], 'r')

codons = {
            "GCT":0, "GCC":1, "GCA":2, "GCG":3,
            "CGT":4, "CGC":5, "CGA":6, "CGG":7, "AGA":8, "AGG":9,
            "AAT":10, "AAC":11,
            "GAT":12, "GAC":13,
            "TGT":14, "TGC":15,
            "CAA":16, "CAG":17,
            "GAA":18, "GAG":19,
            "GGT":20, "GGC":21, "GGA":22, "GGG":23,
            "CAT":24, "CAC":25,
            "ATT":26, "ATC":27, "ATA":28,
            "ATG":29,
            "TTA":30, "TTG":31, "CTT":32, "CTC":33, "CTA":34, "CTG":35,
            "AAA":36, "AAG":37,
            "TTT":38, "TTC":39,
            "CCT":40, "CCC":41, "CCA":42, "CCG":43,
            "TCT":44, "TCC":45, "TCA":46, "TCG":47, "AGT":48, "AGC":49,
            "ACT":50, "ACC":51, "ACA":52, "ACG":53,
            "TGG":54,
            "TAT":55, "TAC":56,
            "GTT":57, "GTC":58, "GTA":59, "GTG":60,
            "TAA":61, "TGA":62, "TAG":63,
            }

fa = {}

for line in seq_file.readlines():
    if line.startswith('>'):
        name = line.strip().replace('>', '')
        fa[name] = ''
    else:
        fa[name] = fa[name] + line.strip()

for ID, seq in fa.items():
    codon_freq = '\t'.join(count_codons(seq, codons))
    gc = seq.replace('C', 'F').replace('G', 'F').count('F')/len(seq)
    print "%s\t%f\t%s" % (ID,gc,codon_freq)







