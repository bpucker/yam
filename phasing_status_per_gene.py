### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###

__usage__ = """
					python phasing_status_per_gene.py
					--cov <COVERAGE_FILE>
					--gff <ANNOTATION_FILE>
					--out <OUTPUT_FOLDER>
					--vcf <VCF_FILE>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import os, glob, sys, gzip, re
import matplotlib.pyplot as plt

# --- end of imports --- #

def get_mean( values ):
	"""! @brief calculate mean of given list of values """
	
	if values != []:
		return sum( values ) / float( len( values ) )
	else:
		return 0


def load_coverage( cov_file ):
	"""! @brief load coverage from given file """
	
	if cov_file[-3:].lower() == "cov" or cov_file[-3:].lower() == "txt":	#uncompressed coverage file
		coverage = {}
		with open( cov_file, "r" ) as f:
			line = f.readline()
			prev_chr = line.split('\t')[0]
			cov = []
			while line:
				parts = line.strip().split('\t')
				if parts[0] != prev_chr:
					coverage.update( { prev_chr: cov } )
					prev_chr = parts[0]
					cov = []
				cov.append( float( parts[2] ) )
				line = f.readline()
			coverage.update( { prev_chr: cov } )
		return coverage
	
	else:	#compressed coverage file
		coverage = {}
		with gzip.open( cov_file, "rb" ) as f:
			line = f.readline()
			prev_chr = line.split('\t')[0]
			cov = []
			while line:
				parts = line.strip().split('\t')
				if parts[0] != prev_chr:
					coverage.update( { prev_chr: cov } )
					prev_chr = parts[0]
					cov = []
				cov.append( float( parts[2] ) )
				line = f.readline()
			coverage.update( { prev_chr: cov } )
		return coverage


def load_genes( gff_file ):
	"""! @brief load gene positions from GFF3 file """
	
	genes = {}
	with open( gff_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[2] == "gene":
					try:
						try:
							try:
								genes.update( { re.findall( "contig\d+\.g\d+", parts[-1] )[0]: { 'chr': parts[0], 'start': int(parts[3]), 'end': int( parts[4] ) } } )
							except IndexError:
								genes.update( { re.findall( "tig\d+\.g\d+", parts[-1] )[0]: { 'chr': parts[0], 'start': int(parts[3]), 'end': int( parts[4] ) } } )
						except IndexError:
							if ';' in parts[-1]:
								genes.update( { parts[-1].split(';')[0].split('=')[1]: { 'chr': parts[0], 'start': int(parts[3]), 'end': int( parts[4] ) } } )
							else:
								genes.update( { parts[-1].split('=')[1]: { 'chr': parts[0], 'start': int(parts[3]), 'end': int( parts[4] ) } } )
					except IndexError:
						print parts[-1]
			line = f.readline()
	return genes


def calculate_cov_per_gene( cov, genes ):
	"""! @brief calculate coverage per gene """
	
	cov_per_gene = {}
	for key in genes.keys():
		values = cov[ genes[ key ]['chr'] ][ genes[ key ]['start']: genes[ key ]['end'] ]
		cov_per_gene.update( { key: get_mean( values ) } )
	return cov_per_gene


def generate_figure( values_to_plot, fig_file, max_value ):
	"""! @brief generate figure with coverage values per gene """
	
	fig, ax = plt.subplots()
	
	ax.hist( values_to_plot, bins=max_value, range=( 0, max_value ) )
	ax.set_xlabel( "average coverage per gene" )
	ax.set_ylabel( "number of genes" )
	
	ax.set_xlim( 0, max_value )
	
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	
	fig.savefig( fig_file, dpi=600 )
	plt.close( "all" )


def load_cov_per_gene( doc_file ):
	"""! @brief load coverage per gene """
	
	coverages = {}
	with open( doc_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			coverages.update( { parts[0]: float( parts[1] ) } )
			line = f.readline()
	return coverages


def load_variants( vcf_file, min_allele_freq, max_allele_freq ):
	"""! @brief load all variants from given VCF file """
	
	variants = {}
	with open( vcf_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if not "," in parts[4]:
					try:
						if parts[-1] != "./.":
							ref, alt = map( int, parts[-1].split(':')[1].split(',') )
							if ( ref+alt ) > 0:
								if min_allele_freq < float( alt ) / ( ref+alt ) < max_allele_freq:
									try:
										try:
											variants[ re.findall( "contig\d+", parts[0] )[0] ].append( int( parts[1] ) )
										except IndexError:
											variants[ re.findall( "tig\d+", parts[0] )[0] ].append( int( parts[1] ) )
									except KeyError:
										try:
											variants.update( { re.findall( "contig\d+", parts[0] )[0]: [ int( parts[1] ) ] } )
										except IndexError:
											variants.update( { re.findall( "tig\d+", parts[0] )[0]: [ int( parts[1] ) ] } )
					except IndexError:
						print parts[-1]
			line = f.readline()
	return variants


def count_variants_per_gene_set( variants, genes ):
	"""! @brief count variants per gene """
	
	variants_per_gene = []
	for key in genes.keys():
		counter = 0
		gene = genes[ key ]
		try:
			for var in variants[ gene['chr'] ]:
				if gene['start'] <= var <= gene['end']:
					counter += 1
		except KeyError:
			pass
		variants_per_gene.append( counter )
	return variants_per_gene


def assess_variant_distribution_in_genes( cov_per_gene, variants, genes, haploid_coverage, filename ):
	"""! @brief assess variant distribution in genes """
	
	# --- classify genes into phased and merged --- #
	phased_genes = {}
	merged_genes = {}
	low_cov_genes = {}
	high_cov_genes = {}
	for gene in cov_per_gene.keys():
		if 0.5*haploid_coverage <= cov_per_gene[ gene ] <= 1.5*haploid_coverage:
			phased_genes.update( {gene: genes[ gene ] } )
		elif 1.5*haploid_coverage < cov_per_gene[ gene ] <= 2.5*haploid_coverage:
			merged_genes.update( {gene: genes[ gene ] } )
		elif cov_per_gene[ gene ] < 0.5*haploid_coverage:
			low_cov_genes.update( {gene: genes[ gene ] } )
		else:
			high_cov_genes.update( {gene: genes[ gene ] } )
	
	print "number of phased genes: " + str( len( phased_genes.keys() ) ) + "\t" + str( 100*len( phased_genes.keys() ) / len( genes.keys() ) ) + "%"
	print "number of merged genes: " + str( len( merged_genes.keys() ) ) + "\t" + str( 100*len( merged_genes.keys() ) / len( genes.keys() ) ) + "%"
	print "number of low coverage genes: " + str( len( low_cov_genes.keys() ) ) + "\t" + str( 100*len( low_cov_genes.keys() ) / len( genes.keys() ) ) + "%"
	print "number of high coverage genes: " + str( len( high_cov_genes.keys() ) ) + "\t" + str( 100*len( high_cov_genes.keys() ) / len( genes.keys() ) ) + "%"
	
	
	# --- count heterozygous variants in phased / merged genes --- #
	variants_per_phased_gene = count_variants_per_gene_set( variants, phased_genes )
	variants_per_merged_gene = count_variants_per_gene_set( variants, merged_genes )
	
	fig, ax = plt.subplots()
	
	ax.boxplot( [ variants_per_phased_gene, variants_per_merged_gene ], labels=[ "phased", "merged" ], showmeans=True )
	ax.set_ylabel( "variants per gene" )
	
	fig.savefig( filename, dpi=600 )
	plt.close( "all" )
	
	print "average number of heterozygous variants in phased genes: " + str( get_mean( variants_per_phased_gene ) )
	print "average number of heterozygous variants in merged genes: " + str( get_mean( variants_per_merged_gene ) )


def main( arguments ):
	"""! @brief runs everything """
	
	cov_file = arguments[ arguments.index( '--cov' )+1 ]
	gff_file = arguments[ arguments.index( '--gff' )+1 ]
	output_dir = arguments[ arguments.index( '--out' )+1 ]
	vcf_file = arguments[ arguments.index( '--vcf' )+1 ]
	
	if '--maxval' in arguments:
		max_value = int( arguments[ arguments.index( '--maxval' )+1 ] )
	else:
		max_value = 150
	
	if '--hapcov' in arguments:
		haploid_coverage = int( arguments[ arguments.index( '--hapcov' )+1 ] )
	else:
		haploid_coverage = 50
	
	
	min_allele_freq = 0.3	#min alt allele frequency to consider variant 'heterozygous'
	max_allele_freq = 0.7	#max alt allele frequency to consider variant 'heterozygous'
	
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	doc_file = output_dir + "value_documentation.txt"
	fig_file = output_dir + "coverage_per_gene.pdf"
	fig_file2 = output_dir + "variants_in_genes.pdf"
	
	genes = load_genes( gff_file )
	if not os.path.isfile( doc_file ):
		cov = load_coverage( cov_file )
		cov_per_gene = calculate_cov_per_gene( cov, genes )
		with open( doc_file, "w" ) as out:
			for gene in sorted( cov_per_gene.keys() ):
				out.write( gene + '\t' + str( cov_per_gene[gene] ) + '\n' )
	else:
		cov_per_gene = load_cov_per_gene( doc_file )
		
	generate_figure( cov_per_gene.values(), fig_file, max_value )
	
	variants = load_variants( vcf_file, min_allele_freq, max_allele_freq )
	assess_variant_distribution_in_genes( cov_per_gene, variants, genes, haploid_coverage, fig_file2 )


if '--cov' in sys.argv and '--gff' in sys.argv and '--out' in sys.argv and '--vcf' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
