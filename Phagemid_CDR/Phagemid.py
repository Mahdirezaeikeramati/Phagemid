from mypackage import translation_Rev3 as TRANS
from mypackage import CDRfinderRev4 as PCDR

####### Main Madule Code
# Modules 
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
import re     #this library is used so that I can use the "search" function
import os     #this is needed for using directory paths and manipulating them 

str =""       #initializing string variable for raw data input
str = input("Enter the path of the data files - please use \ backslash when typing in directory path EXAMPLE==> F:\Bio_Python\codes    :");  #User will enter the name of text file for me
if str!='':
   os.chdir(str) # absolute path, using \ and r prefix
TRANS.Maintrans()
PCDR.CDRfinder()
