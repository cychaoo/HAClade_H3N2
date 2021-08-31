from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
import os,sys

# python mutli-alignment.py "D:\kfcteam-cy\SARS-CoV-2 with Dr.Gong\Flu\gasid_h3n2_seq_processing\raw\gisaid_epiflu_sequence_HA_0813.fasta" "D:\kfcteam-cy\SARS-CoV-2 with Dr.Gong\Flu\gasid_h3n2_seq_processing\raw\metadata_HA_0813.csv"


thisdir = os.path.abspath(os.path.dirname(__file__))
cwd = os.getcwd()

# def main(sysargs = sys.argv[1:]):
#     parser = argparse.ArgumentParser(prog = _program,
#     description='mutli-alignment before imputing gap',
#     usage='''[options]'''

#     parser.add_argument()
#     parser.add_argument()
#     parser.add_argument()

# file with sequences
sequence_file = sys.argv[1]
lineage_file = sys.argv[2]
# f_out = 'output_alignment.fasta'

def substring(seqId, seq, index):
    newSeq = ""
    if index == 0:
        newSeq = newSeq + seq
    else:
        newSeq = newSeq + seq[index:]

    return newSeq


# length must the same
# MAFFT alignment before imputing gap
# then substring "ATG"

def imputeGap(data):
    records = SeqIO.parse(data, 'fasta')
    records = list(records) 
    # substring 'ATG'
    for record in records:
        for i in range(0,len(record),3):
            codon = str(record.seq.upper())[i:i+3]
            if (codon=='ATG'):
                record.seq = Seq(substring(record.id,record.seq,i))
                # print(record.seq)
                break
    # make a copy, otherwise our generator
    # is exhausted after calculating maxlen
    maxlen = max(len(record.seq) for record in records)

    # pad sequences so that they all have the same length
    for record in records:
        if len(record.seq) != maxlen:
            sequence = str(record.seq).ljust(maxlen, '-')
            record.seq = Seq(sequence)
    assert all(len(record.seq) == maxlen for record in records)

    # write to temporary file and do alignment
    output_file = '{}_output.fasta'.format(os.path.splitext(data)[0])
    with open(output_file, 'w') as f:
        SeqIO.write(records, f, 'fasta')
    alignment = AlignIO.read(output_file, "fasta")
    # SeqIO.write(alignment, f_out, 'fasta')

    print(f"Mutli-alignment finshed! \n Output saved:{output_file}")
    return alignment


alignment = imputeGap(sequence_file)
