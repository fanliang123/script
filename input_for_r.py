#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import division
import re
import sys

fasta = open('Sinorhizobium_meliloti_1021.fasta', 'r')
gff = open('Sinorhizobium_meliloti_1021.gff', 'r')
blast_tab = open(sys.argv[1], 'r')
COG_file = open('1021.output.2', 'r')
input_file =  open(sys.argv[1]+'.in', 'w')


fa = {}
cog_fam = {}
cog_cla = {}

for line in fasta.readlines():
    if line.startswith('>'):
        name = line.strip().replace('>', '')
        fa[name] = ''
    else:
        fa[name] = fa[name] + line.strip()

GFF = gff.readlines()

for line in COG_file.readlines():
    if line.startswith('#'):
        continue
    info = line.strip().split("\t")
    cog_fam[info[0]] = info[1]
    cog_cla[info[0]] = info[2]


for line in blast_tab.readlines():
    info = line.strip().split("\t")
    ID = info[0].replace('D_', '')
    for g in GFF:
        if g.startswith('#'):
            continue
        g_info = g.strip().split("\t")
        if g_info[2] != 'CDS':
            continue
        searchObj = re.search(r'.+Name\=(.+?)\;.+product\=(.+?)\;.+', g_info[8])
        # print(ID, searchObj.group(1))
        if ID == searchObj.group(1):
        # if ID == g_info[8].replace('ID=', '').split(';')[0]:
            chrom = g_info[0]
            strand = g_info[6]
            start = g_info[3]
            end = g_info[4]
            seq = fa[chrom][(int(start) - 1) : (int(end))]
            gc = seq.replace('C', 'F').replace('G', 'F').count('F')/len(seq)
            if ID not in cog_fam:
                cog_fam[ID] = 'NA'
                cog_cla[ID] = 'NA'
            info.extend([chrom, strand, start, end, str(gc), cog_fam[ID], cog_cla[ID]])
            input_file.write("\t".join(info)+"\n")
