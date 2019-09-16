import csv
from Bio import SeqIO, Seq, pairwise2

def get_seq(start, strand, fasta='Genome Info/U18466.2.fasta', extend=100):
	for rec in SeqIO.parse(open(fasta), 'fasta'):
		genome = str(rec.seq)
	if strand == '+':
		return genome[start-1-extend : start+extend]
	return str(Seq.reverse_complement(genome[start-1-extend : start+extend]))

def get_overlapped(start, end, strand, tss_poss):
	for tss_pos, tss_strand in tss_poss:
		if tss_strand == strand and start <= tss_pos and end >= tss_pos:
			return tss_pos, tss_strand
	return False

def count_lgaps(seq):
	n = 0
	for x in seq:
		if x == '-':
			n += 1
		else:
			return n
	return n

def get_upstream_aln(seq, ref, start, end, tss, strand, downstream=5):
	alns = pairwise2.align.globalms(seq, ref, 1, -10, -10, -0.1, one_alignment_only=True)
	seq_aln, ref_aln = alns[0][:2]
	if strand == '+':
		ref_start = start - len([x for x in ref_aln.strip('-') if x == '-'])
		seq_start = ref_start
		if seq_aln.startswith('-'):
			seq_start = ref_start + count_lgaps(seq_aln)
		shift = tss - seq_start
	else:
		ref_end = end + len([x for x in ref_aln.strip('-') if x == '-'])
		seq_end = ref_end
		if seq_aln.startswith('-'):
			seq_end = ref_end - count_lgaps(seq_aln)
		shift = seq_end - tss
	
	return seq[:(shift+downstream)][::-1], shift

def count_mismatches(seq, ref, first=25):
	if seq.startswith('G'):
		seq = seq[1:]
		ref = ref[1:]
	return len([i for i in range(first) if seq[i] != ref[i]])

tss_positions = {(int(x.split()[2]), x.split()[-1]): x.split()[3] for x in open('Genome Info/158_pTSS_NA_pTSS_check_NEW_incl_alt_TSSs.bed')}
tss2extendedReads = {}

for file in ['Early-5h/S1_comparison.txt', 'Early-5h/S2_comparison.txt', 
			'Late-16h/S3_comparison.txt', 'Late-16h/S4_comparison.txt']:
	repl = file.split('/')[1].split('_')[0]
	for row in csv.DictReader(open(file), delimiter='\t'):
		start, end, strand = int(row['start']), int(row['end']), row['strand']
		if int(row['difference']) == 0 or count_mismatches(row['sequenced'], row['reference']) == 0:
			continue
		tss = get_overlapped(start, end, strand, tss_positions)
		if tss:
			tss_pos, tss_strand = tss
			if strand == '+':
				shift = tss_pos - start
			else:
				shift = end - tss_pos
			if shift > 1:
				if tss not in tss2extendedReads:
					tss2extendedReads[tss] = []
				tss2extendedReads[tss].append((row['sequenced'], shift))

	for tss in tss2extendedReads:
		tss_pos, tss_strand = tss

		ref = get_seq(tss_pos, tss_strand)[::-1]
		with open('TSS_seq_upstream/'+tss_positions[tss]+'.fa', 'a') as f:
			if repl == 'S1':
				print('>ref', ref, file=f, sep='\n')
			for i, (read, shift) in enumerate(tss2extendedReads[tss]):
				print('>%d_%d_%s'%(i,shift,repl), read[::-1], file=f, sep='\n')





			"""



			#if additional_seq_reversed[-1] == 'G': ## trim the capped G
			#	additional_seq_reversed = additional_seq_reversed[:-1]
			if shift > 1 and shift < 30:
				upstream, real_shift = get_upstream_aln(row['sequenced'], row['reference'], start, end, tss_pos, strand)
				if upstream.endswith('G'):
					upstream = upstream.rstrip('G')
					real_shift -= 1
				if real_shift > 1:
					with open('TSS_seq_upstream/'+tss_positions[tss]+'.fa', 'a') as f:
						print('>%s %d %s %s %s %d'%(repl,tss_pos,start,end,strand,real_shift), 
							row['sequenced'],
							row['reference'],
							upstream, file=f, sep='\n')

			"""

