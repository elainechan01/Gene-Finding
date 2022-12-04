import time
st = time.time()
from exonChainMain import exonChain
geneFile ="./InputSequenceFile/Mycoplasma pneumoniae M129 chromosome.fasta"
genomeFile = "./InputSequenceFile/Guillardia theta nucleomorph chromosome 1, complete sequence.fasta"
geneSeq = exonChain(geneFile, genomeFile)
et = time.time()

print("\n\nExecution time:", et - st)