#!/usr/bin/python

from Genome import whole_genome
from Genome import chrlength

def Nucleosome_code(positionfile, spename):
    # Positionfile is BED format file
    genome = whole_genome(spename)
    (seq_at, seq_cg) = genome.getnucleosomecode()
    at = [0] * 251
    cg = [0] *251
    total = 0
    for line in open(positionfile).xreadlines():
        line = line.strip().split('\t')
	if line[0] not in chrlength[spename].keys():
	    continue
        if line[5] == '+':
            chr = line[0]
            position = int(line[1]) - 1
            begin = position - 50
            if begin < 0:
                continue
            end = position + 201
            if end > chrlength[spename][chr]:
                continue
            total += 1
            subseq_at = seq_at[chr][begin:end]
	    subseq_cg = seq_cg[chr][begin:end]
            for i in range(251):
                at[i] += int(subseq_at[i])
		cg[i] += int(subseq_cg[i])
        elif line[5] == '-':
            chr = line[0]
            position = int(line[2]) - 1
            begin = position - 200
            if begin < 0:
                continue
            end = position + 51
            if end > chrlength[spename][chr]:
                continue
            total += 1
            subseq = ''
            subseqold_at = seq_at[chr][begin:end]
	    subseqold_cg = seq_cg[chr][begin:end]
            for i in reversed(subseqold_at):
                subseq_at += i
	    for i in reversed(subseqold_cg):
		subseq_cg += i
            for i in range(251):
                at[i] += int(subseq_at[i])
		cg[i] += int(subseq_cg[i])
    for i in range(251):
        at[i] = at[i] / float(total)
	cg[i] = cg[i] / float(total)
    return at,cg

if __name__ == '__main__':
    import sys
    from optparse import OptionParser
    usage = 'usage: python Nuc_code_detection.py options positionfile'
    parser = OptionParser(usage)
    
    parser.add_option("-s", "--spename", dest="spename",
                      default="yeast", help="species name" )
    parser.add_option("-o", "--output", dest="output", help="output file")
    (options, args) = parser.parse_args()
    if len(args) < 1:
        parser.print_help()
        sys.exit()
        
    (at,cg) = Nucleosome_code(args[0], options.spename)
    for i in range(251):
        at[i] = str(at[i])
	cg[i] = str(cg[i])
    ofile = open(options.output,'w')
    ofile.write('Code_AT' + ' <- c(' + ', '.join(at) + ')')
    ofile.write('Code_CG' + ' <- c(' + ', '.join(cg) + ')')
    ofile.close()

