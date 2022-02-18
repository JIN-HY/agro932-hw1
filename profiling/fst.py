import pandas as pd
import numpy as np














geno = pd.read_csv("largedata/geno.txt", sep = ",")

def gt2allele(ind):
    ind = [x.split("/")  for x in ind]
    return ind

for i in range(5,407):
    geno.iloc[,i] = gt2allele(geno.iloc[,i])

