## needs bam file with just aligned reads
## and original fastq file

import os
from Bio import SeqIO, pairwise2

r = 'S2'

if '%s_trimmed_mapped_extendedGenome.fasta' % r not in os.listdir('.'):
	cmd = 'samtools view %s_trimmed_mapped.bam | cut -f 1 > %s_trimmed_mapped.names' % (r, r)
	os.system(cmd)
	cmd = 'seqtk subseq %s.fastq.gz %s_trimmed_mapped.names | fastqToFa /dev/stdin %s_trimmed_mapped.fasta' % (r, r, r)
	os.system(cmd)
	cmd = 'bamToBed -i %s_trimmed_mapped.bam > %s_trimmed_mapped.bed' % (r, r)
	os.system(cmd)
	cmd = 'bedtools slop -s -l 25 -i %s_trimmed_mapped.bed -g ../Genome\\ Info/U18466.2.genome_length -r 0 | bedtools getfasta -name+ -fi ../Genome\\ Info/U18466.2.fasta -bed /dev/stdin -s > %s_trimmed_mapped_extendedGenome.fasta' % (r, r)
	os.system(cmd)


read2sequenced = {rec.description.split()[0]: str(rec.seq) for rec in SeqIO.parse(open('%s_trimmed_mapped.fasta' % r), 'fasta')}
read2ref = {rec.description.split('::')[0]: str(rec.seq) for rec in SeqIO.parse(open('%s_trimmed_mapped_extendedGenome.fasta' % r), 'fasta')}
read2pos = {}
for rec in SeqIO.parse(open('%s_trimmed_mapped_extendedGenome.fasta' % r), 'fasta'):
	txt_spl = rec.description.split('::')
	pos = txt_spl[1].split(':')[1]
	start = pos.split('-')[0]
	end = pos.split('-')[1].split('(')[0]
	strand = pos.split('(')[1].split(')')[0]
	read2pos[txt_spl[0]] = (start, end, strand) 

with open('%s_comparison.txt' % r, 'w') as f:
	print('sequenced', 'reference', 'difference', 'start', 'end', 'strand', sep='\t', file=f)
	for read in read2sequenced:
		if read in read2ref:
			seq = read2sequenced[read]
			ref = read2ref[read]
			mismatches = 0
			for i in range(25):
				if seq[i] != ref[i]:
					mismatches += 1
			print(seq, ref, mismatches, *read2pos[read], file=f, sep='\t')

