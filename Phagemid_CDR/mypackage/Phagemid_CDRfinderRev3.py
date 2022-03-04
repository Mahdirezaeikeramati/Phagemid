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

def Phagemid_CDRfinder():
    import string
    n=0
    leftend = 'CAGGTGCAG'
    ### Importer ##########################################################
    output = []
    DNAout = []
    Proteinout = []
    Nameout = []
   # os.chdir(path)
    missing_values = ["n/a", "na", "--"]
    col_list=["Phage ID","NBX","Antigen","Phage ELISA: Antigen A450","Phage ELISA: PBS A450,","Phage Library Lot #","llama #","Lymphocyte Day #","Biopanning  Experiment #","Phage ELISA Experiment #","Current Decision","Protein Sequence","DNA Sequence","Sequencing Data","CDR1","CDR2","CDR3","CDR Family Members for Cloned NBXs","CDR Family ID for Cloned NBXs","Forward Primer","Reverse Primer","Notes","Cloning week","Cloning Initials","Construct ID"]
    data = pd.read_csv("Antibody_Discovery_Table.csv",names=col_list)
    data=data[10:]
    data.columns.name = None
    data1=data[['Phage ID','NBX','Protein Sequence','CDR1','CDR2','CDR3']]
    
#print data1
    data1=data1.dropna()
    df = open("phagemidsequencingProteinout2.txt", "r")
    df = df.read()
    df = df.split('\n')
    del df[-1]
    df = np.array(df)
    df = np.reshape(df,(-1,2))
    Nameout = df[:,[0]]
    Proteinout = df[:,[1]]
    Nameout = np.array(Nameout)
    Proteinout = np.array(Proteinout)
    print("Comparing sequences")
    #Proteinout, indices = np.unique(Proteinout, return_index=True) # remove duplicates the next function can not handle duplicates
    #Nameout = np.take(Nameout, indices) # using the indice form the pervious function to remove duplicates from names
    Nameout = Nameout.flatten()
    Nameout = Nameout.tolist()
    Proteinout = Proteinout.flatten()
    Proteinout = Proteinout.tolist()
    print("Hunting For CDR")
    Misslist = []
    CDRout = []
    CDRTout = []
    CDRTfout=[]
    for x in range (len(Nameout)):
        ProteinIn = str(Proteinout[x])
        NameIn = str(Nameout[x])
        Index = regex.search(r'([TAPSVG][GS][GEDNASQ][SAT][VL][SRTQGKA][LIV][SAFT][CA][ATSVERLQKPG][YAGDASVIPT][SAPTNY]){e<=2}[A-Z]{5,12}([AMLIKVWG][VLSMGADTYRENH][WY][VFYHSINM]R[RQLTEH][RPAGFSTQNVD][TPQAL][GEV][KTNRQESL][APQESRGKDNLV]){e<=2}', ProteinIn)
        if Index ==None:
            VarRegion1 = ''
            Misslist.append(NameIn)
        else:
            VarRegion1 = Index.group()
        VarRegion1 = str(VarRegion1[14:-10])

        Index = regex.search(r'(R[QLTHER][RAVGFSTQNDP][APQTL][GVE][LKNTRQES][GPAESRGKQDNVL][PGRLQV][AVEDIQ][QSFTLWVYIMGR][MLIVA][SAGT]){e<=2}[A-Z]{6,16}(([NYHSFAVR][TVAGSEKNPQYDRL][EDGYQSRPTAN][YFSDTAP][MLVEA][NSVEARPMKIQ][SADVG][GR][FTMSIAW][ISFANTV][LTAIV][SDAWYT][KRSVGNM][EVNYPD][DYVKNGRMHTIAS][GPKDATEVIS][KREDNAVQ][NQTKDSI])){e<=2}', ProteinIn)
        if Index ==None:
            VarRegion2 = ''
            Misslist.append(NameIn)
        else:
            VarRegion2 = Index.group()
        VarRegion2 = str(VarRegion2[15:-19])

        Index = regex.search(r'([YLFSWTDHN][MLV][REQYTLKH][ML][NKSRHDT][TMSKRNYGAD][MLV][MKNSELITRQ][PTLRDSVAF][VGEALD][ED][STA][GDA][GLEIAVD][Y][VNLFTHADSEY][C]){e<=2}[A-Z]{4,24}([KRTQWS][G][QPGEKR][G][ITAS][QEPLHR][V][TINS][VI][ASC][SAQTP]){e<=2}', ProteinIn)
        if Index ==None:
            VarRegion3 = ''
            Misslist.append(NameIn)
        else:
            VarRegion3 = Index.group()
        VarRegion3 = str(VarRegion3[19:-11])

        CDR =  NameIn + ' ' + ProteinIn + ' ' + VarRegion1 + ' ' + VarRegion2 + ' ' + VarRegion3
        #CDR = CDR.split (' ')
        #CDRout.append(CDR)
        STCDR1=start_stop_codon(ProteinIn,VarRegion1)
        STCDR2=start_stop_codon(ProteinIn,VarRegion2)
        STCDR3=start_stop_codon(ProteinIn,VarRegion3)
        Frame1=ProteinIn[0:STCDR1]
        Frame2=ProteinIn[STCDR1+len(VarRegion1):STCDR2]
        Frame3=ProteinIn[STCDR2+len(VarRegion2):STCDR3]
        NBXnew=Frame1 + ' ' +  VarRegion1 + ' ' + Frame2 + ' ' +  VarRegion2 + ' ' + Frame3 + ' ' +  VarRegion3
        CDRT=NameIn + ' ' + Frame1 + ' ' +  VarRegion1 + ' ' + Frame2 + ' ' +  VarRegion2 + ' ' + Frame3 + ' ' +  VarRegion3
        CDRTable1=Frame1+VarRegion1+Frame2 + VarRegion2 + Frame3 +  VarRegion3
        CDRTable=[NameIn , Frame1 ,  VarRegion1 , Frame2 ,  VarRegion2 , Frame3 ,  VarRegion3]
        CDRT = CDRT.split (' ')
        CDRTout.append(CDRT)
        species = ["FR1="+Frame1,"Cdr1="+VarRegion1,"FR2="+Frame2,"Cdr2="+VarRegion2,"FR3="+Frame3,"Cdr3="+VarRegion3]
        count = [len(Frame1),len(VarRegion1),len(Frame2),len(VarRegion2),len(Frame3),len(VarRegion3)]
        NBX1=''
        NBX2=''
        NBX3=''
        FNBX=''
        NBXname=''
        name1=''
        nama=data1["Phage ID"]
        name1==data1.index[nama.loc[1:]== NameIn]
        ## WE are looking for new DNA/RNA Amino Acids Sequences If already it was considered ignor it else continue  just subs if (name1==''): with if (name1!=''):
        if (name1==''):
            NBXname=data1['NBX']
            Cdata=data1['CDR1']
            Bdata=data1['CDR2']
            Ddata=data1['CDR3']
            Fdata=data1['Protein Sequence']
            Maxi=NBXname[max(data1.index)]
            Mamxi=Maxi[6:8]
            NB1=data1.index.values[Cdata.loc[1:]== VarRegion1]
           # print NB1
            NB2=data1.index[Bdata.loc[1:]== VarRegion2]
            NB3=data1.index[Ddata.loc[1:]== VarRegion3]
            FNB=data1.index[Fdata.loc[1:]== CDRTable1]
            if (len(NB1)==0):
                #print ("NB1 is empety")
                X=0
            else:
                NBX1=NB1[0]
                X=1
               # print NBX1
            if (len(NB2)==0):
                #print ("NB1 is empety")
                Y=0
            else:
                NBX2=NB2[0]
                Y=1
                #print NBX2
            if (len(NB3)==0):
                #print ("NB1 is empety")
                Z=0
            else:
                NBX3=NB3[0]
                Z=1
            w=1
            message=''
            Antigen='PVP-'+NameIn[4:7]
            NBXName='' 
           
            if ((X==1) and (Y==1) and (Z==1)):
                if FNB!='':
                   W=1
                   message= ('Ignore '+'duplicated '+'NBX'+' with '+ data1['NBX'][NBX1]+'and sequence'+data1['Phage ID'][NBX3])
                   CDRTf=NameIn + ' ' + NBXName + ' ' + Antigen + ' ' + message+' '
                   CDRTf = CDRTf.split (' ')
                   CDRTfout.append(CDRTf)
                  
                else:
                   w=0
                   if (int(Mamxi)>=100):
                      Ma=""
                   else:
                       Ma=str(int((Maxi[5:6])))
                   n=n+1
                   Alpha=list(string.ascii_lowercase[n-1:n])
                   NBXname= "NBX"+NameIn[1:3]+Ma+str(int(Mamxi)+n)+"_F"+Alpha
                   #print "we find new FR"
            
###  Writing new NBX sequences
            if ((X==0) or (Y==0) or (Z==0)):
                #print "We find New NBX"
                if (int(Mamxi)>=100):
                    Ma=""
                else:
                    Ma=str(int((Maxi[5:6])))
                n=n+1
                NBXName= str("NBX"+NameIn[1:3]+Ma+str(int(Mamxi)+n))
                CDRTf=NameIn + ' ' + NBXName + ' ' + Antigen + ' ' + message + ' ' + Frame1 + ' ' +  VarRegion1 + ' ' + Frame2 + ' ' +  VarRegion2 + ' ' + Frame3 + ' ' +  VarRegion3
                CDRTf = CDRTf.split (' ')
                CDRTfout.append(CDRTf)
            #print "CDR Table Update"
            df = pd.DataFrame(CDRTfout, columns=['Name','NBXName','Antigen','message','FR1','CDR1','FR2','CDR2','FR3','CDR3'])
            df.to_csv("output.csv")
    print "CDR Table was Successfully Updated!"
