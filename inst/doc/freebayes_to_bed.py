import sys


def get_info_map(info):
    info_map = {}
    fields = info.split(';')
    for field in fields:
        if '=' in field:
            idx = field.index('=')
            name = field[:idx]
            value = field[idx+1:]
            info_map[name] = value
    return info_map


def get_counts_info(counts, info_map):
    counts_info = {}
    fields = counts.split(':')
    counts_info['GT'] = fields[info_map.index('GT')]
    counts_info['DP'] = fields[info_map.index('DP')]
    counts_info['RO'] = fields[info_map.index('RO')]
    return counts_info

min_read_depth = float(sys.argv[1])
min_AF = float(sys.argv[2])

pass_count = 0
filter_count = 0

for line in sys.stdin:
    line = line.rstrip()
    if line.startswith('#'):
        continue
    else:
        fields = line.split('\t')
        if len(fields) > 8:
	    info = fields[7]
            info_map = get_info_map(info)
            type = info_map['TYPE']            
            
	    info = fields[8]
            info_map = info.split(':')
	    counts = fields[9]
            counts_info = get_counts_info(counts, info_map)
            
            GT = counts_info['GT']

            DP = 0
            if 'DP' in counts_info:
                DP = abs(float(counts_info['DP']))
            
            RO = 0
            if 'RO' in counts_info:
                RO = abs(float(counts_info['RO']))
            
            RAF = 0
            if DP > 0:
                RAF = RO/DP
            
            MAF = RAF
            if RAF > 0.5:
                MAF = 1 - RAF
            MAF = round(MAF, 3) 
            
            if type == 'snp' and GT == '0/1' and DP > min_read_depth and MAF > min_AF:
                pass_count += 1
                seq = (str(fields[0]).replace('chr', ''), str(int(fields[1])-1), str(int(fields[1])), str(MAF))
                output = '\t'.join(seq)
                print output
            else:
                filter_count += 1

sys.stderr.write('Passed: ' + str(pass_count) + '\n')
sys.stderr.write('Filtered: ' + str(filter_count) + '\n')
sys.stderr.write('Done.\n')

