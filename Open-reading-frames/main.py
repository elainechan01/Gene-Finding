import time
from openReadingFrames import OpenReadingFrame
import pandas as pd

genomeFile = "./datasets/Guillardia theta nucleomorph chromosome 1, complete sequence.fasta"
f = open(genomeFile, "r")
f.readline()
f.readline()
dna = ''
f_contents = f.readlines()
for line in f_contents:
    line = line.rstrip()
    dna += line
dna
st = time.time()
orf = OpenReadingFrame(dna, 300)
orfs = orf.run()
print(orfs)
et = time.time()
print("ORF Execution time: ", et-st)
st = time.time()
orf = OpenReadingFrame(dna, withLikelihoodRatio=True)
orfs, table = orf.run()
print(orfs)
# (pd.DataFrame.from_dict(data=table, orient='index').to_csv('table.csv', header=False))
et = time.time()
print("ORF Execution time w likelihood ratio plot: ", et-st)
f.close()