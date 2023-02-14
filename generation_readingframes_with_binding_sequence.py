#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 16:33:09 2022

@author: friederike
"""

import re
import itertools

#listreadingframes:
minlengthreadingframe=10
maxlengthreadingframe=40
aa_in_core=6
aa_before_binding_site_max = 20
aa_before_binding_site_min = 0
aa_after_binding_site_max = 15
aa_after_binding_site_min = 6
aromatic_amino_acids_forward=('CAT',"CAC","TTT","TTC","TAT","TAC","TGG")
aromatic_amino_acids_reverse=("CCA","GTA","ATA","GAA","AAA","GTG","ATG")
filenamereadingframes="Output/readingframes_only_common_binding_sites.txt"

startcodon1 = ('ATG', "GTG")
stopcodon1=('TGA', 'TAA', 'TAG')

startcodon2 = ("CAT", "CAC")
stopcodon2=("TCA","TTA","CTA","GTA")


#source=https://parts.igem.org/Ribosome_Binding_Sites/Catalog ,http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html
ribosomalbindingsites1=["AAAGA[A-Z]{3}GA[A-Z]{3}","AACTAGAATCACCTCTTGCTTTTGGGTAAGAAAGAGGAGA","AACTAGAATCACCTCTTGGATTTGGGTATTAAAGAGGAGA","ATTAAAGAGGAGAAA","TTCACACAGGAAACC",'GTGTG','GTGTGTCTAG','TCACACAGGAAACCGGTTCGATG','TCACACAGGAAAGGCCTCGATG','TCACACAGGACGGCCGGATG','TCTCACGTGTGTCAAG','TCTCACGTGTGT','CATCCCT','TCACATCCCT','TCACATCCCTCC','ACTGCACGAGGTAACACAAG','TACGAGGAGGATGAAGAGTA','ACTTTACTTATGAGGGAGTA','ACGAAGACGGAGACTTCTAA','AACCCTCAGGAGGTAAACCA','AAGACATGGAGACACATTTA','ACTGCACGAGGTAACACAAG','GAGGAGGATGAAGAGTAATGTGAAGAGCTG','AGGAGGTCATC','GCAAGCTCTTTTTTCAGTTGTCTC','CTGATAGTTAAAATCACCAGCATGA','TAAAAACAAGAGGAAAACAA','TCTCCTCTTT','ACGGAGAAGCAGCGAA','GAGGTTGGGACAAG','TAAATGTATCCGTTTATAAGGACAGCCCGA','CTCTTAAGTGGGAGCGGCT','CTCTACCGGAGAAATT','CTCATCGTTAAAGAGCGACTAC','CTCAGCCTGTACCTGGAGAGCCTTTC','CTCAAGGAGG','GAGAGG','AGGAGGATTACAA','AAAGAGGAGAAA','TCACACAGGAAAG','GGAAGAGG','TTTCTCCTCTTTAAT','TCACACAGGAAAGGCCTCG','ATTAAAGAGGAGAAATTAAGC',' TCGTTTCTGAAAAATTTTCGTTTCTGAAAA','TGGCTAACATAGGGT','TGGCTAACTGAGGAT','TGGCTAACCCAGGGT','TGGCTAACTCAGGTG','TGGCTAACCCTGGTA','TGGCTAACTTGGGAC','TGGCTAACGCAGGTC','TGGCTAACATCGGTG','TTAATTAAGGAAAAGATCT','CAGAAGAGGATATTAATA','TTGATAAGGAATTGTA','TCAGAGGAGATAATTTA','TGACACGTTGAGCGGTATGA','ACAGATAACAGGAGTAAGTA','TAAAGGGAGAAAAAT','GAGTCTTGAGGTAACTAT','TCAGGAATATTAAAAACGCT','ATTTGAAGGAAAATATT','CAAAAACATACTGCAGGAAT','TGCCATTGCAAAGGAGAAGACT','AAGGGGGAATTCAAAT','AAGGGGTGCAGAAT','AGGTGGAATCACAG','ATAGATAAAAATGGTAACAAT','GGGATATAGCCTGAGGGGCCTGTA','CGGCAATAACAGAGGCGATTT','ATTAAAGAGGAGAAATA','TCACACAGGAAAGTA','AAAGGAGGTGT','AGAGGTGGTGT','AGGAGG','GAGG','TAAAGGAGGAA','AAAGGTGGTGAA','AGGAAACAGAACC','ATATTAAGAGGAGGAG','AGAGAACAAGGAGGGG','GATTGGGATAAATAAT','ATCAACCGGGGTACAT','TTTGGAGATTTTCAAC','AAAAAAGGTAATTCAA','CATAAGGTAATTCACA','ATAAGGAGTCTTAATC','GTTCCGGCTAAGTAAC','TAATGGAAACTTCCTC','TCGCTGGGGGTCAAAG','ATTTGAGGGGGATTCA','AATTTAGGTCAGAAG','AATCAATAGGAGAAATCAAT','TTAAAGAGGAGAAATACTAG']
ribosomalbindingsites2=["[A-Z]{3}TC[A-Z]{3}TCTTT","TCTCCTCTTTCTTACCCAAAAGCAAGAGGTGATTCTAGTT","TCTCCTCTTTAATACCCAAATCCAAGAGGTGATTCTAGTT","TTTCTCCTCTTTAAT",'GGTTTCCTGTGTGAA','CACAC','CTAGACACAC','CATCGAACCGGTTTCCTGTGTGA','CATCGAGGCCTTTCCTGTGTGA','CATCCGGCCGTCCTGTGTGA','CTTGACACACGTGAGA','ACACACGTGAGA','AGGGATG','AGGGATGTGA','GGAGGGATGTGA','CTTGTGTTACCTCGTGCAGT','TACTCTTCATCCTCCTCGTA','TACTCCCTCATAAGTAAAGT','TTAGAAGTCTCCGTCTTCGT','TGGTTTACCTCCTGAGGGTT','TAAATGTGTCTCCATGTCTT','CTTGTGTTACCTCGTGCAGT','CAGCTCTTCACATTACTCTTCATCCTCCTC','GATGACCTCCT','GAGACAACTGAAAAAAGAGCTTGC','TCATGCTGGTGATTTTAACTATCAG','TTGTTTTCCTCTTGTTTTTA','AAAGAGGAGA','TTCGCTGCTTCTCCGT','CTTGTCCCAACCTC','TCGGGCTGTCCTTATAAACGGATACATTTA','AGCCGCTCCCACTTAAGAG','AATTTCTCCGGTAGAG','GTAGTCGCTCTTTAACGATGAG','GAAAGGCTCTCCAGGTACAGGCTGAG','CCTCCTTGAG','CCTCTC','TTGTAATCCTCCT','TTTCTCCTCTTT','CTTTCCTGTGTGA','CCTCTTCC','ATTAAAGAGGAGAAA','CGAGGCCTTTCCTGTGTGA','GCTTAATTTCTCCTCTTTAAT','TTTTCAGAAACGAAAATTTTTCAGAAACGA','ACCCTATGTTAGCCA','ATCCTCAGTTAGCCA','ACCCTGGGTTAGCCA','CACCTGAGTTAGCCA','TACCAGGGTTAGCCA','GTCCCAAGTTAGCCA','GACCTGCGTTAGCCA','CACCGATGTTAGCCA','AGATCTTTTCCTTAATTAA','TATTAATATCCTCTTCTG','TACAATTCCTTATCAA','TAAATTATCTCCTCTGA','TCATACCGCTCAACGTGTCA','TACTTACTCCTGTTATCTGT','ATTTTTCTCCCTTTA','ATAGTTACCTCAAGACTC','AGCGTTTTTAATATTCCTGA','AATATTTTCCTTCAAAT','ATTCCTGCAGTATGTTTTTG','AGTCTTCTCCTTTGCAATGGCA','ATTTGAATTCCCCCTT','ATTCTGCACCCCTT','CTGTGATTCCACCT','ATTGTTACCATTTTTATCTAT','TACAGGCCCCTCAGGCTATATCCC','AAATCGCCTCTGTTATTGCCG','TATTTCTCCTCTTTAAT','TACTTTCCTGTGTGA','ACACCTCCTTT','ACACCACCTCT','CCTCCT','CCTC','TTCCTCCTTTA','TTCACCACCTTT','GGTTCTGTTTCCT','CTCCTCCTCTTAATAT','CCCCTCCTTGTTCTCT','ATTATTTATCCCAATC','ATGTACCCCGGTTGAT','GTTGAAAATCTCCAAA','TTGAATTACCTTTTTT','TGTGAATTACCTTATG','GATTAAGACTCCTTAT','GTTACTTAGCCGGAAC','GAGGAAGTTTCCATTA','CTTTGACCCCCAGCGA','TGAATCCCCCTCAAAT','CTTCTGACCTAAATT','ATTGATTTCTCCTATTGATT','CTAGTATTTCTCCTCTTTAA']

#create different shine-dalgarno sites:
dalgarno1=[]
dalgarno2=[]

#also look for common cores
commoncore=["TCCTGGTACATCTGGTACTGA","GCCTGGTACCTCTGGTACTGA","GCCTGGTACATCTGGTACTGA","GCCTGGTACCACTGGTACTAA"]
commoncorerev=["TCAGTACCAGATGTACCAGGA","TCAGTACCAGAGGTACCAGGC","TCAGTACCAGATGTACCAGGC","TTAGTACCAGTGGTACCAGGC"]
common_binding_site = ["AAAAGCCTCAAG", "AAAAGCCTGAAG", "AAAAGCCTGAAG", "AAAAGCTTGAAG", "AAAAGTCTGAAG", "AAGAGCCGTAAG", "AAGAGCCTCAAA", "AAGAGCCTCAAG", "AAGAGCCTGAAA", "AAGAGCCTGAAG", "AAGAGCCTGAAG", "AAGAGCTTGAAG", "AAGAGTCTCAAA", "AAGAGTCTCAAG", "AAGAGTCTGAAA", "AAGAGTCTGAAG",
                       "AAGCCGCTGAAG", "AAGCCGTTGAAG", "AAGTCCCTGAAG", "AGGTCTCGGAAG", "CCCAGCCTCAAG", "CCCAGCCTGAAG", "CCCAGTCTGAAG", "CCGAGCCTCAAG", "CCGAGCCTGAAG", "CCGAGCCTGAAG", "CGCTCCCTCAAG", "CGCTCGCTCAAG", "CGGAGCCTGAAG", "CGGAGCCTGAAG", "CGTTCGCTCAAG", "CGTTCGCTGAAA", "GAAAGCCTCAAG"]
common_binding_site_reversed = ["CTTACGGCTCTT", "CTTCAACGGCTT", "CTTCAAGCTCTT", "CTTCAAGCTTTT", "CTTCAGACTCTT", "CTTCAGACTGGG", "CTTCAGACTTTT", "CTTCAGCGGCTT", "CTTCAGGCTCCG", "CTTCAGGCTCCG", "CTTCAGGCTCGG", "CTTCAGGCTCGG", "CTTCAGGCTCTT", "CTTCAGGCTCTT", "CTTCAGGCTGGG", "CTTCAGGCTTTT",
                               "CTTCAGGCTTTT", "CTTCAGGGACTT", "CTTCCGAGACCT", "CTTGAGACTCTT", "CTTGAGCGAACG", "CTTGAGCGAGCG", "CTTGAGGCTCGG", "CTTGAGGCTCTT", "CTTGAGGCTGGG", "CTTGAGGCTTTC", "CTTGAGGCTTTT", "CTTGAGGGAGCG", "TTTCAGACTCTT", "TTTCAGCGAACG", "TTTCAGGCTCTT", "TTTGAGACTCTT", "TTTGAGGCTCTT"]
def generate_necessary_motifs():
    listreadingframe = []
    forward_aa_list=list(itertools.product(aromatic_amino_acids_forward, repeat=2))
    forward_aa=[]
    for pair in forward_aa_list:
        forward_aa+=["".join(pair)]
    reverse_aa_list=itertools.product(aromatic_amino_acids_reverse, repeat=2)
    reverse_aa=[]
    for pair in reverse_aa_list:
        reverse_aa+=["".join(pair)]
    for i in range(minlengthreadingframe,maxlengthreadingframe):
        i=i*3
        for end in stopcodon1:
            for n in range (1,5):
                for start in startcodon1:
                    for forward in forward_aa:
                        listreadingframe=listreadingframe+[start+"[A-Z]{"+str(i)+"}"+"[A-Z]{"+str((4-n)*3)+"}"+str(forward)+"[A-Z]{"+str((n)*3)+"}"+str(end)]
    for i in range(minlengthreadingframe,maxlengthreadingframe):
        i=i*3
        for end in stopcodon2:
            for n in range (1,5):
                for start in startcodon2:
                    for reverse in reverse_aa:
                        listreadingframe=listreadingframe+[end+"[A-Z]{"+str((4-n)*3)+"}"+str(reverse)+"[A-Z]{"+str((n)*3)+"}"+"[A-Z]{"+str(i)+"}"+str(start)]
    return listreadingframe
#listreadingframe=generate_necessary_motifs()
listreadingframe = []
#create reading frame patterns

for aa_after_binding_site in range(aa_after_binding_site_min, aa_after_binding_site_max+1):
    bases_after_binding_site=aa_after_binding_site*3
    for end in stopcodon1:
        for aa_before_binding_site in range (aa_before_binding_site_min, aa_after_binding_site_max +1):
            bases_before_binding_site = aa_before_binding_site *3
            for start in startcodon1:
                for p450_binding_site in common_binding_site:
                    listreadingframe = listreadingframe + \
                        [start+"[A-Z]{"+str(bases_before_binding_site)+"}" + p450_binding_site +
                         "[A-Z]{"+str(bases_after_binding_site)+"}"+str(end)]
for aa_after_binding_site in range(aa_after_binding_site_min, aa_after_binding_site_max+1):
    bases_after_binding_site = aa_after_binding_site*3
    for end in stopcodon2:
        for aa_before_binding_site in range(aa_before_binding_site_min, aa_after_binding_site_max + 1):
            bases_before_binding_site = aa_before_binding_site * 3
            for start in startcodon2:
                for p450_binding_site in common_binding_site_reversed:
                    listreadingframe = listreadingframe + \
                        [end+"[A-Z]{"+str(bases_after_binding_site)+"}" + p450_binding_site +
                         "[A-Z]{"+str(bases_before_binding_site)+"}"+str(start)]
# for aa_after_binding_site in range(minlengthreadingframe-7,maxlengthreadingframe):
#     aa_after_binding_site=aa_after_binding_site*3
#     for end in commoncore:
#         for start in startcodon1:
#             listreadingframe=listreadingframe+[start+"[A-Z]{"+str(aa_after_binding_site)+"}"+str(end)]
# for aa_after_binding_site in range(minlengthreadingframe-7,maxlengthreadingframe):
#     aa_after_binding_site=aa_after_binding_site*3
#     for end in commoncorerev:
#         for start in startcodon2:
#             newframe=end+"[A-Z]{"+str(aa_after_binding_site)+"}"+start
#             listreadingframe=listreadingframe+ [newframe]
#safe file
with open(filenamereadingframes, 'w') as f:
    for s in listreadingframe:
        f.write(str(s) + '\n')