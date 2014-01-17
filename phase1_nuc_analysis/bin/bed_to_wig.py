#!/usr/bin/python

import sys

chrlength = {}
chrlength['yeast'] = {'chr1' : 230208,  'chr2' : 813178,  'chr3' : 316617,  'chr4' : 1531918,
                          'chr5' : 576869,  'chr6' : 270148,  'chr7' : 1090947, 'chr8' : 562643,
                          'chr9' : 439885,  'chr10' : 745745,  'chr11' : 666454,  'chr12' : 1078175,
                          'chr13' : 924429,  'chr14' : 784333,  'chr15' : 1091289, 'chr16' : 948062,
                          'chrM' : 85779}

chrlength['dm3'] = {'chr2L':23011544,'chr2LHet':368872,'chr2R':21146708,'chr2RHet':3288761,
                    'chr3L':24543557,'chr3LHet':2555491,'chr3R':27905053,'chr3RHet':2517507,
		    'chr4':1351857,'chrX':22422827,'chrXHet':204112,'chrYHet':347038,'chrU':10049037,
		    'chrUextra':29004656,'chrM':19517}

chrlength['zv9'] = { 'chr1':60348388,'chr2':60300536,'chr3':63268876,'chr4':62094675,
'chr5':75682077,'chr6':59938731,'chr7':77276063,'chr8':56184765,'chr9':58232459,'chr10':46591166,'chr11':46661319,
'chr12':50697278,'chr13':54093808,'chr14':53733891,'chr15':47442429,'chr16':58780683,'chr17':53984731,'chr18':49877488,
'chr19':50254551,'chr20':55952140,'chr21':44544065,'chr22':42261000,'chr23':46386876,'chr24':43947580,'chr25':38499472}

chrlength['hg19'] = {'chr1':249250621,'chr2':243199373,'chr3':198022430,'chr4':191154276,'chr5':180915260,'chr6':171115067,
'chr7':159138663,'chrX':155270560,'chr8':146364022,'chr9':141213431,'chr10':135534747,'chr11':135006516,'chr12':133851895,
'chr13':115169878,'chr14':107349540,'chr15':102531392,'chr16':90354753,'chr17':81195210,'chr18':78077248,'chr20':63025520,
'chrY':59373566,'chr19':59128983,'chr22':51304566,'chr21':48129895}


def tag2hash(bedfile, spename, extension, shift):
    values = {}
    for chrname in chrlength[spename].keys():
        values[chrname] = [0] * chrlength[spename][chrname]
        
    for line in open(bedfile).xreadlines():
        line = line.strip().split()
        if line[0] not in chrlength[spename].keys():
            continue
            
        if line[5] == '+':           # Tag in plus strand
            b = int(line[1]) + shift
            e = b + extension
            for k in xrange(max(b, 1), min(e, chrlength[spename][line[0]] + 1)):
                values[line[0]][k - 1] += 1
                            
        elif line[5] == '-':        # Tag in minus strand
            e = int(line[2]) - shift
            b = e - extension
            for k in xrange(max(b, 1), min(e, chrlength[spename][line[0]] + 1)):
                values[line[0]][k - 1] += 1
        else:
            continue
    return values

def tag_to_wig(hash_values, wigfile, step, spename):
    
    namelist = chrlength[spename].keys()
    namelist.sort()
    for name in namelist:
        fileout = open(name + '_' + wigfile + '.wig','w')
	print >>fileout, "track type=wiggle_0 name=" + name
        print >>fileout, "fixedStep  chrom=" + name + "  start=1  step=" + str(step)
        for k in xrange(0, len(hash_values[name]), step):
            print >>fileout, hash_values[name][k]
    fileout.close()
    
if __name__ == '__main__':
    import sys
    from optparse import OptionParser
    
    usage = "usage: python bed_to_wig.py -i input bed file -o output wig file"
    parser = OptionParser(usage)
    
    parser.add_option("-i",dest="input",help="input bed file")
    parser.add_option("-o",dest="output",help="output wig file")
    parser.add_option("-s",dest="species",type="string",default='zv9',help="choose the species")
    parser.add_option("--step",dest="step",type="int",default=10,help="step of output wig")
    parser.add_option("--extension",dest="extension",type="int",default=73,help="The extension of bed reads")
    parser.add_option("--shift",dest="shift",type="int",default=36,help="The shift of bed reads")
    
    (options,args) = parser.parse_args()
    
    if not options.input or not options.output:
        parser.print_help()
        sys.exit()
    
    values = tag2hash(options.input, options.species, options.extension, options.shift)
    tag_to_wig(values,options.output,options.step,options.species)
    
    
    
    
    
