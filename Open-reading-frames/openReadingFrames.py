from collections import Counter
import pandas as pd
from copy import deepcopy

class OpenReadingFrame:
    """Class to implement Open Reading Frames"""
    def __init__(self, dnaSequence: str, threshold: int = 0, withLikelihoodRatio: bool=False) -> None:
        """Constructor Method
        
        Required Args
            dnaSequence (str): dna sequence
        
        Optional Args
            threshold (int): threshold value to determine candidate ORFs
            withLikelihoodRatio (bool): implement likelihood ratio plot statistical approach
        """
        self.codon_table =  {
            'TTT': 'F',     'CTT': 'L',     'ATT': 'I',     'GTT': 'V',
            'TTC': 'F',     'CTC': 'L',     'ATC': 'I',     'GTC': 'V',
            'TTA': 'L',     'CTA': 'L',     'ATA': 'I',     'GTA': 'V',
            'TTG': 'L',     'CTG': 'L',     'ATG': 'M',     'GTG': 'V',
            'TCT': 'S',     'CCT': 'P',     'ACT': 'T',     'GCT': 'A',
            'TCC': 'S',     'CCC': 'P',     'ACC': 'T',     'GCC': 'A',
            'TCA': 'S',     'CCA': 'P',     'ACA': 'T',     'GCA': 'A',
            'TCG': 'S',     'CCG': 'P',     'ACG': 'T',     'GCG': 'A',
            'TAT': 'Y',     'CAT': 'H',     'AAT': 'N',     'GAT': 'D',
            'TAC': 'Y',     'CAC': 'H',     'AAC': 'N',     'GAC': 'D',
            'TAA': 'STOP',  'CAA': 'Q',     'AAA': 'K',     'GAA': 'E',
            'TAG': 'STOP',  'CAG': 'Q',     'AAG': 'K',     'GAG': 'E',
            'TGT': 'C',     'CGT': 'R',     'AGT': 'S',     'GGT': 'G',
            'TGC': 'C',     'CGC': 'R',     'AGC': 'S',     'GGC': 'G',
            'TGA': 'STOP',  'CGA': 'R',     'AGA': 'R',     'GGA': 'G',
            'TGG': 'W',     'CGG': 'R',     'AGG': 'R',     'GGG': 'G'
        }
        self.dna = dnaSequence
        self.withLikelihoodRatio = withLikelihoodRatio
        self.codonFrequencyTable = {}
        self.nucleotideFrequencyTable = {}
        self.nucleotidePositionLikelihoodPlot = {}
        self.threshold = threshold

    def initCodonFrequencyTable(self) -> None:
        """Method to initialize nucleotide frequency table by calculating the probability of the occurence of the nucleotide in the dna sequence"""
        occurrenceFromStart = Counter([self.dna[i:i+3] for i in range(len(self.dna)-2)])
        dna = self.findComplement()
        occurrenceFromEnd = Counter([dna[i:i+3] for i in range(len(dna)-2)])
        occurrenceFull = dict(list(occurrenceFromStart.items()) + list(occurrenceFromEnd.items())) 
        # determine the occurence count of nucleotides of the dna string
        for codon in self.codon_table.keys():
            self.codonFrequencyTable[codon] = 1       # laplaces rule
        for item in occurrenceFull:
            self.codonFrequencyTable[item] += occurrenceFull[item]
        self.codonUsage = deepcopy(self.codonFrequencyTable)
        for codon, occurence in self.codonFrequencyTable.items():
            self.codonFrequencyTable[codon] = occurence/sum(self.codonFrequencyTable.values())
            self.codonUsage[codon] = [occurence, occurence/sum(self.codonFrequencyTable.values())]

    def translateCodon(self, codon: str) -> str:
        """Method to translate codons to protein

        Required Args
            codon (str): codon of base triplets

        Returns
            str: protein based on codon_table
        """
        try:
            return self.codon_table[codon]
        except KeyError:
            return ''

    def findComplement(self) -> str:
        """Method to translate reverse string based on nucleotides' complements

        Returns
            str: reverse dna sequence
        """
        complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
        return ''.join([complement[nucleotide] for nucleotide in reversed(self.dna)])

    def findAllORF(self, reverse: bool=False, returnRawCodon: bool=False) -> list:
        """Method to find all Open Reading Frames that starts with a Start codon, ATG and ends with a Stop codon, TAA, TAG or TGA

        Optional Args
            reverse (bool): find ORF of reverse dna sequence
            returnRawCodon (bool): return non-translated codons instead of translated concatenated ORFs

        Returns
            list: codons (translated and concatenated if returnRawCodon is specified)
        """
        if reverse:
            frames = {
                -1: {},
                -2: {},
                -3: {}
            }
        else:
            frames = {
                1: {},
                2: {},
                3: {}
            }
        # determine type of string to evaluate
        if reverse:
            dna = self.findComplement()
        else:
            dna = self.dna

        possibleORFs = []
        
        # find all start codons ('ATG') and store its position
        start_codons = {}
        # frame 1
        for i in range(len(dna)):
            if self.translateCodon(dna[i:i+3]) == 'M':
                start_codons[i] = dna[i:i+3]
                if reverse:
                    frames[-1][i] = 0
                else:
                    frames[1][i] = 0
        # frame 2
        dna = self.dna[1:]
        for i in range(len(dna)-1):
            if self.translateCodon(dna[i:i+3]) == 'M':
                start_codons[i] = dna[i:i+3]
                if reverse:
                    frames[-2][i] = 0
                else:
                    frames[2][i] = 0
        # frame 2
        dna = self.dna[2:]
        for i in range(len(dna)-2):
            if self.translateCodon(dna[i:i+3]) == 'M':
                start_codons[i] = dna[i:i+3]
                if reverse:
                    frames[-3][i] = 0
                else:
                    frames[3][i] = 0
        
        # find relative stop codons for all start codons ('TAA', 'TAG', 'TGA')
        for start in start_codons.keys():
            for index in range(start+3, len(dna), 3):
                if self.translateCodon(dna[index:index+3]) == 'STOP':
                    for frame in frames:
                        if start in frames[frame].keys():
                            frames[frame][start] = index+3
                    # add ORF by splicing the dna sequence starting from the position of the start codon and ending with the position of the stop codon
                    if not returnRawCodon:
                        # return the translated and concatenated version of ORFs if returnRawCodon is not specified
                        possibleORFs.append(''.join([self.translateCodon(dna[i:i+3]) for i in range(start,index,3)]))
                    else:
                        # return only the raw codons if returnRawCodon is specified
                        possibleORFs.append([dna[i:i+3] for i in range(start,index,3)])
                    break
        return possibleORFs, frames

    def findMostLikelyORFs(self, sequences: list) -> list:
        """Method to find all Open Reading Frames that starts with a Start codon, ATG and ends with a Stop codon, TAA, TAG or TGA

        Required Args
            sequences (list): list of ORFs as raw codons

        Returns
            list: codons (translated and concatenated if returnRawCodon is specified)
        """
        # initialize score table and result array of most likely ORFs
        score = {}
        # determine likelihood ratio for each ORF
        # the score is computed as the log-likelihood ratio of the probability of the codon at its position over the frequency of occurence in the genome (likelihood of occurence at current position over other random positions in the dna sequence)
        for seq in sequences:
            score[''.join([self.translateCodon(codon) for codon in seq])] = sum([self.codonFrequencyTable[codon] for codon in seq])
        # return most probable ORF
        for seq in score:
            if score[seq] == max(score.values()):
                return seq, score[seq]

    def run(self) -> list:
        """Run Method

        Returns
            list: translated ORFs
        """
        # condition: likelihood ratio plot is preferred, return all ORFs that are most likely to occur
        if self.withLikelihoodRatio:
            self.initCodonFrequencyTable()
            possibleSequencesFromStart, framesStart = self.findAllORF(returnRawCodon=True)
            possibleSequencesFromEnd, framesEnd = self.findAllORF(reverse=True, returnRawCodon=True)
            mostLikelySeq, score = self.findMostLikelyORFs(possibleSequencesFromStart + possibleSequencesFromEnd)
            return mostLikelySeq, score, self.codonUsage
        # condition: normal ORF is preferred, return all ORFs that exceed a threshold value
        else:
            possibleSequencesFromStart, framesStart = self.findAllORF()
            possibleSequencesFromEnd, framesEnd = self.findAllORF(reverse=True)
            possibleSequencesFull = {seq: len(seq) for seq in set(possibleSequencesFromStart + possibleSequencesFromEnd) if len(seq) >= self.threshold}
            for seq in possibleSequencesFull:
                if possibleSequencesFull[seq] == max(possibleSequencesFull.values()):
                    return seq, len(seq)
