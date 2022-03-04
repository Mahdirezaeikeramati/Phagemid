
import regex
import numpy as np
import pandas as pd
from Bio import SeqIO, pairwise2
from Bio.pairwise2 import format_alignment  #optional
from termcolor import colored
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from bokeh.plotting import figure, show, output_file #outils graphiques
## Translate DNA/mRNA to Amino Acid  
def translate(seq): 
    tabletest = [ "TTT" , "CTT", "ATT" , "GTT" ,
           "TTC" , "CTC" , "ATC" , "GTC" ,"TTA" , "CTA" , "ATA" , "GTA" ,
           "TTG", "CTG" , "ATG", "GTG" ,"TCT" , "CCT", "ACT", "GCT",
           "TCC", "CCC", "ACC", "GCC","TCA", "CCA", "ACA", "GCA",
           "TCG", "CCG", "ACG", "GCG" ,"TAT", "CAT", "AAT", "GAT" ,
           "TAC" , "CAC" , "AAC" , "GAC","TAA","CAA", "AAA", "GAA" ,
           "TAG" , "CAG" , "AAG", "GAG","TGT", "CGT", "AGT", "GGT" ,
           "TGC", "CGC" , "AGC" , "GGC","TGA", "CGA", "AGA", "GGA",
           "TGG", "CGG", "AGG", "GGG"]
   # print len(tabletest)
    table = { "TTT" : "F", "CTT" : "L", "ATT" : "I", "GTT" : "V",
           "TTC" : "F", "CTC" : "L", "ATC" : "I", "GTC" : "V",
           "TTA" : "L", "CTA" : "L", "ATA" : "I", "GTA" : "V",
           "TTG" : "L", "CTG" : "L", "ATG" : "M", "GTG" : "V",
           "TCT" : "S", "CCT" : "P", "ACT" : "T", "GCT" : "A",
           "TCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
           "TCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
           "TCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
           "TAT" : "Y", "CAT" : "H", "AAT" : "N", "GAT" : "D",
           "TAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
           "TAA" : "STOP", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
           "TAG" : "STOP", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
           "TGT" : "C", "CGT" : "R", "AGT" : "S", "GGT" : "G",
           "TGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
           "TGA" : "STOP", "CGA" : "R", "AGA" : "R", "GGA" : "G",
           "TGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G"   
           }
    protein1 =""
    
    if len(seq)%3 == 0:
        X=True
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            if codon in tabletest:
#                print "Hi"
##            for j in range(63):
##                 if (codon ==tabletest[j]):
##                     X=True
##                     #print k
##            if X==True:
                 protein1+= table[codon]
            else:
                 X=False
        if X==True:
            protein=protein1
        else:
            protein=""        
    return protein 

def read_seq(inputfile): 
    with open(inputfile, "r") as f: 
      seq = f.read()
      seq = seq.replace("\n", "")
      seq = seq.replace("\r", "") 
    return seq

#find Start and Stop Codon
def start_stop_codon(record,querystring):
    ##print dna
    aln = pairwise2.align.localms(record, querystring,
                                   5, -4,      # match, mismatch
                                   -2, -0.5)   # open, extend

  #  [print(format_alignment(*a)) for a in aln]     #optional

    coords = [i for i in range(len(aln[0][0])) if aln[0][0][i] == aln[0][1][i]]
    subseq = aln[0][0][coords[0]:coords[-1]+1]
    seqstart = str(record).index(subseq)

    #print('Subseq coordinates: {} -> {}'.format(seqstart, seqstart+len(subseq)))
    return seqstart
