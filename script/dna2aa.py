import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio import AlignIO
from Bio import pairwise2
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
import pandas as pd
import os,sys

# file with sequences
sequence_file = sys.argv[1]
lineage_file = sys.argv[2]

def getDataLine(seqId, seq):
    dataLine = []
    dataLine.append(seqId)

    newSeq = ""
    newSeq = newSeq + seq

    dataLine.append(newSeq)

    return dataLine


# altered DNA to protein (DNA to AA)
def translate(seq,metadata):
    metadata = pd.read_csv(metadata,encoding='ISO-8859-1')
    seqs = SeqIO.parse(seq, 'fasta')
    dataList = []
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
        '---': '-',}
    protein =""
    for record in seqs:
        for line in range(len(record)):
            if record.id == metadata['strain'][line]:
                for i in range(0, len(record.seq), 3):
                    if str(record.seq.upper())[i:i + 3] in table:
                        codon = str(record.seq.upper())[i:i + 3]
                        protein+= table[codon]
                    else:
                        codon = "X"
                        protein+= codon
            dataList.append(getDataLine(record.id,Seq(protein)))
            break
    return dataList

tran_alignment = translate(sequence_file,lineage_file)

# export transform data
pd.DataFrame(tran_alignment, columns=['id','Seq']).to_csv('{}_trans_alignment.csv'.format(os.path.splitext(sequence_file)[0]), index = False,encoding='utf-8')

path = os.path.splitext(sequence_file)[0]
print(f'Transforming finshed! \n Output saved: {path}_tran_alignment.csv')