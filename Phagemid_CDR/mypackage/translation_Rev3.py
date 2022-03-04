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
from Bio import SeqIO, pairwise2
from Bio.pairwise2 import format_alignment  #optional
import glob
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
def start_codon(record,querystring):
    ##print dna
    aln = pairwise2.align.localms(record, querystring,
                                   5, -4,      # match, mismatch
                                   -2, -0.5)   # open, extend

  #  [print(format_alignment(*a)) for a in aln]     #optional

    coords = [i for i in range(len(aln[0][0])) if aln[0][0][i] == aln[0][1][i]]
    subseq = aln[0][0][coords[0]:coords[-1]+1]
 #   print(coords)
 #   print(subseq)
    if (subseq!="CAGGTGCAG"):
       seqstart=0 
    else:
        seqstart = str(record).index(subseq)
 #       print(seqstart) 
    #print('Subseq coordinates: {} -> {}'.format(seqstart, seqstart+len(subseq)))
    return seqstart
def stop_codon(record,querystring):
    ##print dna
    aln = pairwise2.align.localms(record, querystring,
                                   5, -4,      # match, mismatch
                                   -2, -0.5)   # open, extend

  #  [print(format_alignment(*a)) for a in aln]     #optional

    coords = [i for i in range(len(aln[0][0])) if aln[0][0][i] == aln[0][1][i]]
    subseq = aln[0][0][coords[0]:coords[-1]+1]
#    print(coords)
#    print(subseq)
    if (subseq!="ACCGTCTCC"):
       seqstart=0 
    else:
        seqstart = str(record).index(subseq)
#        print(seqstart) 
    #print('Subseq coordinates: {} -> {}'.format(seqstart, seqstart+len(subseq)))
    return seqstart

def is_subsequence(lst1, lst2):
    """
        *   Finds if a list is a subsequence of another.

        *   Params:
            *   `lst1` (`list`): The candidate subsequence.
            *   `lst2` (`list`): The parent list.
        *   Return:
            *   (`bool`): A boolean variable indicating whether `lst1` is a subsequence of `lst2`.
    """
    l1, l2 = len(lst1), len(lst2)
    if l1 > l2:  #`l1` must be <= `l2` for `lst1` to be a subsequence of `lst2`.
        return False
    i = j = 0
    d1, d2 = l1, l2
    while i < l1 and j < l2:
        while lst1[i] != lst2[j]:
            j += 1
            d2 -= 1
            if d1 > d2:  #At this point, `lst1` cannot a subsequence of `lst2`.
                return False
        i, j, d1, d2 = i+1, j+1, d1-1, d2-1
        if d1 > d2:
            return False
    return True
## Main code
def Maintrans():
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
    from Bio import SeqIO, pairwise2
    from Bio.pairwise2 import format_alignment  #optional
    import glob
    print ("Reading Sequence Files")
    path=(glob.glob("*.ab1"))
    f = open("phagemidsequencingProteinout2.txt", "w")
    for i in range(len(path)):
          handle = open(path[i], "rb")
          records = list(SeqIO.parse(handle, "abi"))
#print("Found %i records" % len(records))
##print("The last record")
          last_record= records[-1]
          Name=str(path[i][0:14])
     # print(Name)
      #print(repr(last_record.seq))
##print(len(last_record))
          first_record= records[0]
          dna = first_record.seq
     # print dna
          sub1='AGGTGCAG'   #"CAGGTGCAGCT"  
          record=dna
          sub2="ACCGTCTCC"    #"TCACCGTCTCC"
          sub3="TAA"
          sub4="TAG"
          sub5="TGA" 
          mm=start_codon(dna,sub1)
          nn1=stop_codon(dna,sub2)
##          nn2=start_stop_codon(dna,sub3)
##          nn3=start_stop_codon(dna,sub4)
##          nn4=start_stop_codon(dna,sub5)
##      if (nn1<=mm):
##          nn1=10000
##      if (nn2<=mm):
##          nn2=10000
##      if (nn3<=mm):
##          nn3=10000
##      if (nn4<=mm):
##          nn4=10000   
##      nn=min(nn1,nn2,nn3,nn4)
     # print nn1,nn2,nn3,nn4,mm
      
      #print mm,nn+len(sub2)
          Re=(nn1+len(sub2)- mm) % 3
      #print Re
          if Re == 0:
         
             p1 = translate(dna[mm:nn1+len(sub2)])
             print(p1)
             Index=is_subsequence('STOP',p1)
             print(Index)
         #Index = regex.search(r'([STOP])', p1)
         #print Index
             if Index==False:
                 p=p1
          #   print 'hello'
             else:
                 p=""
          else:
             p=""
          if p!="":
             f.write(">")
             f.write(Name)
             f.write("\n")
             f.write(p)
             print(p)
             f.write("\n")
    print(Index)         
    print("Detection of Codon")
    print("Translate Sequence to Amino Acid Sequence") 
    f.close()
