import os

gff = 'U18466.2.gff3'
out = 'genes.bb'

id2name = {}
for line in open(gff):
	tab = line.strip().split('\t')
	if len(tab) < 8:
		continue
	desc = tab[8]
	desc_spl = desc.split(';')
	if desc.startswith('ID') and 'Name=' in desc:
		id2name[desc.split('ID=')[1].split(';')[0]] = desc.split('Name=')[1].split(';')[0]

colors = ['31,119,180', '44,160,44', '31,119,180', '255,127,14', '214,39,40']
feature2col = {}

with open(out, 'w') as f:
	for line in open(gff):
		tab = line.strip().split('\t')
		if len(tab) < 8:
			continue
		desc = tab[8]
		extra = ''
		if 'Note=' in desc:
			extra = desc.split('Note=')[1].split(';')[0].replace('%3B', ';').replace('%2C', ',')[:255]
		if tab[2] not in ['gene', 'region']:
			if 'Parent=' in desc:
				xid = id2name[desc.split('Parent=')[1].split(';')[0]]
			elif 'Name=' in desc:
				xid = desc.split('Name=')[1].split(';')[0]
			else:
				xid = desc.split('ID=')[1].split(';')[0]
			feature = tab[2]
			if feature not in feature2col:
				feature2col[feature] = colors.pop()
			col = feature2col[feature]
			print(tab[0], tab[3], tab[4], xid, 1000, tab[6], tab[3], tab[4], col, extra, sep='\t', file=f)

os.system('sort -k1,1 -k2,2n %s > tmp2' % out)
#os.system('echo "track name=Genes type=bedDetail" > tmp1')
#os.system('cat tmp1 tmp2 > tmp3')
os.system('bedToBigBed -as=genes_table.as -extraIndex=name,geneProduct -type=bed9+1 -tab tmp2 chrom.sizes %s' % out)
#os.system('rm tmp2')
