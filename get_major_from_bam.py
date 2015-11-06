#! /usr/bin/python
import sys, pysam


sam=sys.argv[1]
samfile=pysam.Samfile(sam,'rb')


def select_read(read):
    if read.alignment.mapq>=20 and \
       read.alignment.flag in [147,99,163,83,1107,1187,1123,1171]:
        return read
    else:
        pass

def get_major(bases):
    major ={0:'A',1:'C',2:'G',3:'T'}
    counts=[bases.count('A'),
            bases.count('C'),
            bases.count('G'),
            bases.count('T')]
    return major[counts.index(max(counts))] 


def printseq(seq,output):
    """Print out sequence"""
    print >> output, '>chrM'
    for i in range(0,len(seq),70):
        print >> output, ''.join(seq[i:i+70])


sequence=['N']*16569
outfasta=open(sam+'.major.fa','w+')

for pileupcolumn in samfile.pileup( 'chrM', 0, 16569, stepper='all', max_depth=10000000, mask=False):
    pos_in_ref=pileupcolumn.pos
    bases_in_position=[]
    for pileupread in pileupcolumn.pileups:
        if select_read(pileupread):
            alignment=pileupread.alignment
            pos_in_read=pileupread.qpos
            if ord(alignment.qual[pos_in_read])-33>=30:
               bases_in_position.append(alignment.seq[pos_in_read])
    sequence[pos_in_ref]=get_major(bases_in_position)

printseq(sequence,outfasta)
