import pandas as pd
import numpy as np






geno = pd.read_csv("largedata/geno.txt", sep = ",")

def gt2allele(ind):
    ind = [x.split("/")  for x in ind]
    return ind

def sumalt(ind_call):
    if ind_call[0] == ".": continue
    n_alt = ind_call

for i in range(5,407):
    geno.iloc[,i] = gt2allele(geno.iloc[,i])

