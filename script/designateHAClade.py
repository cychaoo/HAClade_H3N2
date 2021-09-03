import pandas as pd
import os,sys
from Bio import SeqIO
import csv

# flie
transed = sys.argv[1] # amino acids trans alignment
orginal = sys.argv[2] # nucleotide data

with open(transed, newline='') as f:
    reader = csv.reader(f)
    dna = list(reader)
    
dna.pop(0) # remove column name


def getDataLine(seqId, seq):
    dataLine = []
    dataLine.append(seqId)

    newSeq = ""
    newSeq = newSeq + str(seq)

    dataLine.append(newSeq)

    return dataLine


nucList = []
seq_dict = {rec.id : rec.seq.upper() for rec in SeqIO.parse(orginal, "fasta")}
for key in seq_dict.keys():
    nucList.append(getDataLine(key, seq_dict[key]))

newData = pd.DataFrame(dna).rename(columns={0: 'id', 1:'seq_HA'},inplace = False).merge(pd.DataFrame(nucList).rename(columns={0: 'id', 1:'seq_nuc'},inplace = False),how = 'inner' , on = 'id')


def Clade(ID,HAline,nucline):
    # HA1 position +16，但從0開始所以要-1，[15:16] or [15]
    # HA2 position +343
    newLine = []
    newLine.append(ID)
    clade = ""

    if (str(HAline[160:161])=='S') and (str(HAline[174:175])=='F') and (str(HAline[75:76])=='K') and (str(HAline[214:215])=='S') and (str(HAline[238:239])=='I') and (str(HAline[327:328])=='S') and (str(HAline[500:501])=='N') and (str(nucline[1194:1195])=='C') and (str(nucline[1670:1671])=='G'):
        clade = '3b' 
    # 3b
    # HA1 145 S, HA1 159 F, HA1 60 K, 
    # HA1 198 S, HA1 223 I, HA1 312 S, 
    # HA2 158 N, nuc 1195 C, nuc 1671 G

    elif (str(HAline[63:64])=='I') and (str(HAline[60:61])=="N") and (str(nucline[455:456])=='T'):
        clade = '3c'
    # 3c
    # HA1 48 I, HA1 45 N, nuc 456 T

    elif (str(HAline[502:503])=='N') and (str(nucline[692:693])=='A') and (str(nucline[1517:1518])=='G'):
        clade = '3c2'
    # 3c2
    # HA2 160 N, nuc 693 A, nuc 1518 G

    elif (str(HAline[143:144])=='A') and (str(HAline[157:158])=='G') and (str(nucline[1295:1296])=='A'): 
        clade = '3c3'
    # 3c3 
    # HA1 128 A, HA1 142 G, nuc 1296 A

    elif (str(HAline[174:175])=='Y') and (str(HAline[18:19])=='I') and (str(HAline[159:160])=='S') and (str(HAline[175:176])=='T') and (str(HAline[174:175])=='S') and (str(HAline[240:241])=='D') and (str(HAline[153:154])=='S') and (str(HAline[341:342])=='R') and (str(nucline[1259:1260])=='A'):
        clade = '3c2.A'
    #3c2.A	
    # HA1 159 Y, HA1 3 I, HA1 144 S, 
    # HA1 160 T, HA1 159 S, HA1 225 D, 
    # HA1 138 S, HA1 326 R, nuc 1260 A

    elif (str(HAline[98:99])=='R') and (str(HAline[276:277])=='Q') and (str(HAline[77:78])=='K') and (str(HAline[360:361])=='K'):
        clade = '3c3.B'
    # 3c3.B	
    # HA1 83 R, HA1 261 Q, HA1 62 K

    elif (str(HAline[186:187])=='K') and (str(HAline[419:420])=='V') and (str(HAline[497:498])=='E'):
        clade = 'A1'
    # A1 
    # HA1 171 K, HA2 77 V, HA2 155 E

    elif (str(HAline[186:187])=='K') and (str(HAline[419:420])=='V') and (str(HAline[492:493])=='E') and (str(HAline[497:498])=='E') and (str(nucline[80:81])=='A') and (str(nucline[113:114])=='T') and (str(nucline[1483:1484])=='A'):
        clade = 'A1a'
    # A1a	
    # HA1 171 K, A1a HA2 77 V, HA2 150 E,
    # HA2 155 E, nuc 81 A, nuc 114 T, nuc 1484 A

    elif (str(HAline[107:108])=='R') and (str(HAline[326:327])=='Q') and (str(nucline[263:264])=='G') and (str(nucline[537:538])=='C'):
        clade = 'A1b'
    # A1b
    # HA1 92 R, HA1 311 Q, nuc 264 G, nuc 538 C

    elif (str(nucline[327:328])=='A'):
        clade = 'A1b/94N'
    # A1b/94N
    # nuc 328 A

    elif (str(HAline[146:147])=='K') and (str(HAline[77:78])=='G') and (str(HAline[157:158])=='G'):
        clade = 'A1b/131K'
    # A1b/131K	
    # HA1 131 K, HA1 62 G, HA1 142 G

    elif (str(HAline[150:151])=='K') and (str(HAline[77:78])=='G') and (str(HAline[157:158])=='G'):
        clade = 'A1b/135K'
    # A1b/135K
    # HA1 135 K, HA1 62 G, HA1 142 G

    elif (str(HAline[150:151])=='N') and (str(nucline[80:81])=='G'):
        clade = 'A1b/135N'
    # A1b/135N	
    # HA1 135 N, nuc 81 G

    elif (str(HAline[150:151])=='K') and (str(HAline[77:78])=='G') and (str(HAline[157:158])=='G') and (str(HAline[208:209])=='S'):
        clade = 'A1b/137F'
    #A1b/137F	
    # HA1 135 K, HA1 62 G, HA1 142 G,
    #HA1 193 S

    elif (str(HAline[150:151])=='K') and (str(HAline[201:202])=='D') and (str(HAline[205:206])=='N'):
        clade = 'A1b/186D'
    # A1b/186D
    # HA1 135 K, HA1 186 D, HA1 190 N

    elif (str(HAline[146:147])=='K') and (istr(HAline[77:78])=='G') and (str(HAline[157:158])=='G') and (str(HAline[212:213])=='R'):
        clade = 'A1b/197R'
    # A1b/197R	
    # HA1 131 K, HA1 62 G, HA1 142 G,
    # HA1 197 R

    elif (str(HAline[174:175])=='N') and (str(HAline[175:176])=='I'):
        clade = 'A1b/159N'
    # A1b/159N	
    # HA1 159 N, HA1 160 I

    elif (str(HAline[201:202])=='S') and (str(HAline[213:214])=='P'):
        clade = 'A1b/186S'
    # A1b/186S	
    # HA1 186 S, HA1 198 P

    elif (str(HAline[276:277])=='Q') and (str(HAline[157:158])=='K') and (str(nucline[1484:1485])=='T'):
        clade = 'A2'
    # A2	
    # HA1 261 Q, HA1 142 K, nuc 1485 T

    elif (str(nucline[1688:1689])=='T') and (str(nucline[1124:1125])=='A'):
        clade = 'A2/re'   
    # A2/re	
    # nuc 1689 T, nuc 1125 A

    elif (str(HAline[136:137])=='K') and (str(nucline[1133:1134])=='G') and (str(nucline[1319:1320])=='T'):
        clade = 'A3'
    # A3	
    # HA1 121 K, nuc 1134 G, nuc 1320 T

    elif (istr(HAline[207:208])=='T') and (str(HAline[212:213])=='H') and (str(HAline[46:47])=='S') and (str(HAline[68:69])=='N') and (str(HAline[159:160])=='R'):
        clade = 'A4'
    # A4
    # HA1 192 T, HA1 197 H, HA1 31 S,
    # HA1 53 N, HA1 144 R

    else:
        clade = 'NA'

    newLine.append(clade)
    return newLine



# find columns in the data list which clade designation
def findHAClade(data):

    #nucdata = SeqIO.parse(nuc, 'fasta')
    newList = []

    for i in data.index:
        newLine = Clade(data["id"][i],data['seq_HA'][i],data['seq_nuc'][i])
        newList.append(newLine)

    return newList


newDataList = findHAClade(newData)


pd.DataFrame(newDataList).rename(columns={0: 'id', 1:'clade'},inplace = False).to_csv('{}_clade_designation.csv'.format(os.path.splitext(transed)[0]), index = False,encoding='utf-8')

path = os.path.splitext(transed)[0]
print(f'HAClade designation finshed! \n Output saved: {path}_clade_designation.csv')