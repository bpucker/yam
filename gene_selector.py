### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

import sys, re, os
from operator import itemgetter

# --- end of imports --- #

__usage__ = """ python gene_selector.py
							--list <HOMOLOGY_BASED_ANALYSIS_RESULTS (TXT)>
							--gff <INPUT_FILE (FASTA)>
							--out <OUTPUT_FOLDER>
							
							optional:
							--fasta <MATCHING_ASSEMBLY_FILE>
							--name <NAME_OF_OUTPUT_FILES>["x"]
							--sr <FLOAT, score ratio cutoff>[0.1]
							--cov <FLOAT, percentage of query aligned>[0.1]
							
							feature requests and bug reports: bpucker@cebitec.uni-bielefeld.de
						"""

def load_genes_of_interest( list_input_file, sr_cutoff, cov_cutoff ):
	"""! @brief load interesting genes from given file """
	
	genes = {}
	with open( list_input_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if float( parts[1] ) >= sr_cutoff and float( parts[2] ) >= cov_cutoff:
				ID = re.findall( "tig\d+\.g\d+", parts[0] )[0]
				try:
					genes[ ID ]
				except KeyError:
					genes.update( { ID: None } )
			line = f.readline()
	return genes


def load_annotation_per_gene( gff_file ):
	"""! @brief load annotation per gene """
	
	anno_per_gene = {}
	gene_order = []
	with open( gff_file, "r" ) as f:
		line = f.readline()
		while line:
			content = []
			if line[ : len( "# start gene " ) ] == "# start gene ":
				ID = re.findall( "tig\d+\.g\d+", line )[0]
				while line[ : len( "# end gene " ) ] != "# end gene ":
					content.append( line )
					line = f.readline()
				content.append( line )
				anno_per_gene.update( { ID: content } )
				gene_order.append( ID )
			line = f.readline()
	return anno_per_gene, gene_order


def main( arguments ):
	"""! @brief run everything """
	
	list_input_file = arguments[ arguments.index('--list')+1 ]
	gff_file = arguments[ arguments.index('--gff')+1 ]
	output_directory = arguments[ arguments.index('--out')+1 ]
	
	if '--fasta' in arguments:
		fasta_file = arguments[ arguments.index('--fasta')+1 ]
	else:
		fasta_file = ""
	
	if '--name' in arguments:
		name = arguments[ arguments.index('--name')+1 ]
	else:
		name = "x"
	
	if '--sr' in arguments:
		sr_cutoff = float( arguments[ arguments.index('--sr')+1 ] )
	else:
		sr_cutoff = 0.1	#score ratio db/self
	
	if '--cov' in arguments:
		cov_cutoff = float( arguments[ arguments.index('--cov')+1 ] )
	else:
		cov_cutoff = 0.1	#fraction of query covered by alignment
	
	
	perl_script = "/vol/cluster-data/bpucker/bin/scripts/getAnnoFasta.pl"
	
	if not os.path.exists( output_directory ):
		os.makedirs( output_directory )
	
	# --- select genes of interest based on comparison to databases --- #
	genes_of_interest = load_genes_of_interest( list_input_file, sr_cutoff, cov_cutoff )
	print "number of selected genes: " + str( len( genes_of_interest.keys() ) )
	
	# --- load annotation per gene and keep gene order --- #
	annotation_per_gene, gene_order = load_annotation_per_gene( gff_file )
	
	# --- generate filtered GFF3 file --- #
	new_gff_file = output_directory + name + ".gff"
	with open( new_gff_file, "w" ) as out:
		out.write( "##gff-version 3\n#\n" )
		for gene in gene_order:
			try:
				genes_of_interest[ gene ]
				out.write( "".join( annotation_per_gene[ gene ] ) + "#\n" )
			except KeyError:
				pass
	
	# --- extract annotated sequences --- #
	if len( fasta_file ) > 0:
		os.chdir( output_directory )
		cmd = " ".join( [ perl_script, "--seqfile=" + fasta_file, new_gff_file ] )
		os.popen( cmd )


if '--list' in sys.argv and '--gff' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
