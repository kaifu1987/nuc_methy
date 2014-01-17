#!/usr/bin/python

# Yong Zhang @ DFCI

import sys, time
#from numpy import *

chrlength = {}

chrlength['ecoli'] = {'U00096' : 4639675}
chrlength['yeast_UCSC'] = {'chr1' : 230208,  'chr2' : 813136,  'chr3' : 316613,  'chr4' : 1531914,
                      'chr5' : 576869,  'chr6' : 270148,  'chr7' : 1090944, 'chr8' : 562639,
                      'chr9' : 439885,  'chr10' : 745446,  'chr11' : 666445,  'chr12' : 1078173,
                      'chr13' : 924430,  'chr14' : 784328,  'chr15' : 1091285, 'chr16' : 948060,
                      'chrM' : 85779}
chrlength['yeast'] = {'chrI' : 230208,  'chrII' : 813178,  'chrIII' : 316617,  'chrIV' : 1531918,
'chrV' : 576869,  'chrVI' : 270148,  'chrVII' : 1090947, 'chrVIII' : 562643,
'chrIX' : 439885,  'chrX' : 745745,  'chrXI' : 666454, 'chrXII' : 1078175,
'chrXIII' : 924429,  'chrXIV' : 784333,  'chrXV' : 1091289, 'chrXVI' : 948062,
'chrM' : 85779} 
chrlength['yeast_SGD'] = {'chr01' : 230208,  'chr02' : 813178,  'chr03' : 316617,  'chr04' : 1531918,
                      'chr05' : 576869,  'chr06' : 270148,  'chr07' : 1090947, 'chr08' : 562643,
                      'chr09' : 439885,  'chr10' : 745745,  'chr11' : 666454,  'chr12' : 1078175,
                      'chr13' : 924429,  'chr14' : 784333,  'chr15' : 1091289, 'chr16' : 948062,
                      'chrmt' : 85779}
chrlength['zv9'] = { 'chr1':60348388,'chr2':60300536,'chr3':63268876,'chr4':62094675,
'chr5':75682077,'chr6':59938731,'chr7':77276063,'chr8':56184765,'chr9':58232459,'chr10':46591166,'chr11':46661319,
'chr12':50697278,'chr13':54093808,'chr14':53733891,'chr15':47442429,'chr16':58780683,'chr17':53984731,'chr18':49877488,
'chr19':50254551,'chr20':55952140,'chr21':44544065,'chr22':42261000,'chr23':46386876,'chr24':43947580,'chr25':38499472}

import time
time.asctime()
class Genome(object):
    
    def __init__(self, filename):
        temp = []
        for line in open(filename).xreadlines():
            if line[0] == '>':
                continue
            temp.append(line.strip())
        self.seq = ''.join(temp)
    def length(self):
        return len(self.seq)
    
    def ATcontent(self, length, width = 146):
        at = []
        for k in range(length - width):
            subseq = self.seq[k : k + width]
            numat = subseq.count('A') + subseq.count('T') + subseq.count('a') + subseq.count('t')
            numcg = subseq.count('C') + subseq.count('G') + subseq.count('c') + subseq.count('g')
            if (numat + numcg) >= 0.8 * width:
                at.append(str(k + width / 2) + '\t' + str(int(1000 * numat / (numat + numcg)) / 10.0))
        return at

lastfix = {'ecoli' : '.fna', 'yeast' : '.fa', 'yeast_SGD' : '.fsa', 'yeast_UCSC' : '.fa','zv9':'.fa'}

dirname = '/u/home/k/kaifu/project-mcdb/nuc_methy/genome_infor/chromFa'

class whole_genome(object):
    def __init__(self, spename):
        if not chrlength.has_key(spename):
            print 'Wrong species name:'
            sys.exit()
        self.seq = {}
        for chrname in chrlength[spename].keys():
            temp = []
            for line in open(dirname + '/' + chrname + lastfix[spename]).xreadlines():
                if line[0] == '>':
                    continue
                temp.append(line.strip().upper())
            self.seq[chrname] = ''.join(temp)
            
    def getseq(self):
        return self.seq
    def getnucleosomecode(self):
        code_at = {}
        code_cg = {}
        for chrname in self.seq.keys():
	    code_at[chrname] = []
            code_cg[chrname] = []
	    for i in range(len(self.seq[chrname]) - 1):
                pa = self.seq[chrname][i:i+2]
                if pa == 'AT' or pa == 'AA' or pa == 'TA' or pa == 'TT':
                    code_at[chrname].append('1')
                else:
                    code_at[chrname].append('0')
		if pa == 'CG' or pa == 'CC' or pa == 'GC' or pa == 'GG':
		    code_cg[chrname].append('1')
		else:
		    code_cg[chrname].append('0')
	    code_at[chrname].append('0')
	    code_cg[chrname].append('0')

        return code_at, code_cg

def test():
#    dirname = '/home/zy/collaboration/Nuc_Struhl/Data/Genome/ecoli/'
#    fo = open('/home/zy/collaboration/Nuc_Struhl/Analysis/ATcontent/Ecoli_AT.wig', 'w')
#    for name in chrlength['ecoli'].keys():
#        genome = Genome(dirname + name + '.fna')
#        at = genome.ATcontent(chrlength['ecoli'][name])
#        print >>fo, "track type=wiggle_0"
#        print >>fo, "variableStep  chrom=" + name
#        for k in xrange(len(at)):
#            print >>fo, at[k]
#    fo.close()
    
    dirname = '/home/zy/collaboration/Nuc_Struhl/Data/Genome/yeast_SGD/'
    fo = open('/home/zy/collaboration/Nuc_Struhl/Analysis/ATcontent/Yeast_SGD_AT_1k.wig', 'w')
    for name in chrlength['yeast_SGD'].keys():
        genome = Genome(dirname + name + '.fsa')
        length = genome.length()
        print name + ' : ' + str(length)
        at = genome.ATcontent(length, width= 1000)
        print >>fo, "track type=wiggle_0"
        print >>fo, "variableStep  chrom=" + name
        for k in xrange(len(at)):
            print >>fo, at[k]
    fo.close()

def test_convergent_genes():
    gene = Gene('/Users/yzhang/Bioinfo/Nuc_Struhl/Genes_features/Annotation/SGD_gene_SGD.bed')
    gene.convergent_genes('/Users/yzhang/Bioinfo/Nuc_Struhl/Genes_features/Annotation/SGD_gene_SGD_5end_convergent.bed', end = 5)
    gene.convergent_genes('/Users/yzhang/Bioinfo/Nuc_Struhl/Genes_features/Annotation/SGD_gene_SGD_3end_convergent.bed', end = 3)
    gene.share_ends_genes('/Users/yzhang/Bioinfo/Nuc_Struhl/Genes_features/Annotation/SGD_gene_SGD_5end_shared_others_5end.bed', end = 5, end_others = 5)
    gene.share_ends_genes('/Users/yzhang/Bioinfo/Nuc_Struhl/Genes_features/Annotation/SGD_gene_SGD_5end_shared_others_3end.bed', end = 5, end_others = 3)
    gene.share_ends_genes('/Users/yzhang/Bioinfo/Nuc_Struhl/Genes_features/Annotation/SGD_gene_SGD_3end_shared_others_5end.bed', end = 3, end_others = 5)
    gene.share_ends_genes('/Users/yzhang/Bioinfo/Nuc_Struhl/Genes_features/Annotation/SGD_gene_SGD_3end_shared_others_3end.bed', end = 3, end_others = 3)
    
    gene = Gene('/Users/yzhang/Bioinfo/Nuc_Struhl/Genes_features/Annotation/SGD_gene_UCSC.bed')
    gene.convergent_genes('/Users/yzhang/Bioinfo/Nuc_Struhl/Genes_features/Annotation/SGD_gene_UCSC_5end_convergent.bed', end = 5)
    gene.convergent_genes('/Users/yzhang/Bioinfo/Nuc_Struhl/Genes_features/Annotation/SGD_gene_UCSC_3end_convergent.bed', end = 3)
    gene.share_ends_genes('/Users/yzhang/Bioinfo/Nuc_Struhl/Genes_features/Annotation/SGD_gene_UCSC_5end_shared_others_5end.bed', end = 5, end_others = 5)
    gene.share_ends_genes('/Users/yzhang/Bioinfo/Nuc_Struhl/Genes_features/Annotation/SGD_gene_UCSC_5end_shared_others_3end.bed', end = 5, end_others = 3)
    gene.share_ends_genes('/Users/yzhang/Bioinfo/Nuc_Struhl/Genes_features/Annotation/SGD_gene_UCSC_3end_shared_others_5end.bed', end = 3, end_others = 5)
    gene.share_ends_genes('/Users/yzhang/Bioinfo/Nuc_Struhl/Genes_features/Annotation/SGD_gene_UCSC_3end_shared_others_3end.bed', end = 3, end_others = 3)
    
if __name__ == '__main__':
    #test_convergent_genes()
    test()
    
