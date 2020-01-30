### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
				python find_single_copy_phased_orthologs.py
				--phasing <PHASING_FILE>
				--orthogroups <ORTHOGROUP_FILE>
				
				bug reports and feature reqests: bpucker@cebitec.uni-bielefeld.de
				"""

import re, sys, os

# --- end of imports --- #


def load_phased_genes( phasing_status_file, min_cov, max_cov ):
	"""! @brief identify phased genes based on their average coverage """
	
	phased_genes = {}
	with open( phasing_status_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if min_cov < float( parts[1] ) < max_cov:
				phased_genes.update( { parts[0]: None } )
			line = f.readline()
	return phased_genes


def count_phased_single_copy_orthologs( phased_genes, orthogroup_file ):
	"""! @brief count phased single copy orthologs """
	
	counter = 0
	one_to_two = 0
	with open( orthogroup_file, "r" ) as f:
		f.readline()
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			IDs = re.findall( "contig\d+\.g\d+", parts[2] )
			if len( IDs ) == 2:
				one_to_two += 1
				try:
					phased_genes[ IDs[0] ]
					phased_genes[ IDs[1] ]
					counter += 1
				except KeyError:
					pass
			line = f.readline()

	print "1:2 relationships: " + str( one_to_two )
	print "number of phased single copy orthologs: " + str( counter )
	return counter


def main( arguments ):
	"""! @brief run everything """
	
	phasing_status_file = arguments[ arguments.index( '--phasing' )+1 ]
	orthogroup_file = arguments[ arguments.index( '--orthogroups' )+1 ]
	
	min_cov = 75
	max_cov = 150
	
	phased_genes = load_phased_genes( phasing_status_file, min_cov, max_cov )
	
	print "number of phased genes: " + str( len( phased_genes ) )
	
	phased_single_copy_orthologs = count_phased_single_copy_orthologs( phased_genes, orthogroup_file )


if '--phasing' in sys.argv and '--orthogroups' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
