### Boas Pucker ###
### v0.1 ###
### bpucker@cebitec.uni-bielefeld.de ###

__usage__ = """ python match_alleles.py
						--prefix <FULL_PATH_TO_DIRECTORY_FOR_TMP_DATA_AND_RESULTS>
						--cds <FULL_PATH_TO_INPUT_FILE1>
				
						feature requests and bug reports: bpucker@cebitec.uni-bielefeld.de
					"""

import re, os, sys
from operator import itemgetter


# --- end of imports --- #

def load_results_from_BLAST_result_file( BLAST_result_file, cutoff=0.99 ):
	"""! @brief load data from BLAST result file """
	
	data = {}
	with open( BLAST_result_file, "r" ) as f:
		line = f.readline()
		prev_query = line.split('\t')[0]
		hits = []
		while line:
			parts = line.strip().split('\t')
			if parts[0] != prev_query:
				sorted_hits = sorted( hits, key=itemgetter( 'score' ) )
				if len( sorted_hits ) > 2:
					if ( sorted_hits[-3]['score'] / sorted_hits[-2]['score'] ) < cutoff:
						data.update( { sorted_hits[-2]['query']: sorted_hits[-2]['subject'] } )
				elif len( sorted_hits ) == 2:
					data.update( { sorted_hits[-2]['query']: sorted_hits[-2]['subject'] } )
				else:
					pass
				hits = []
				prev_query = parts[0]
			hits.append( { 'query': parts[0], 'subject': parts[1], 'score': float( parts[-1] ) } )
			line = f.readline()
		sorted_hits = sorted( hits, key=itemgetter( 'score' ) )
		if len( sorted_hits ) > 2:
			if ( sorted_hits[-3]['score'] / sorted_hits[-2]['score'] ) > cutoff:
				data.update( { sorted_hits[-2]['query']: sorted_hits[-2]['subject'] } )
		elif len( sorted_hits ) == 2:
			data.update( { sorted_hits[-2]['query']: sorted_hits[-2]['subject'] } )
		else:
			pass
	print "entries in data: " + str( len( data.keys() ) )
	return data


def load_multiple_fasta_file( fasta_file ):
	"""!@brief load content of multiple fasta file """
	
	content = {}
	with open( fasta_file, "r" ) as f:
		header = f.readline().strip()[1:]
		line = f.readline()
		seq = ""
		while line:
			if line[0] == '>':
				content.update( { header: seq } )
				header = line.strip()[1:]
				seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		content.update( { header: seq } )
	return content


def identify_RBHs( data, RBH_file ):
	"""! @brief identify RBHs between alleles """
	
	black_list = {}
	for key in data.keys():
		try:
			value = data[ key ]
			if data[ value ] == key:
				try:
					black_list[ key ]
				except KeyError:
					black_list.update( { value: None } )
		except KeyError:
			pass
	print "number of matched alleles: " + str( len( black_list.keys() ) )


def main( parameters ):
	"""! @brief identifies RBHs between given data sets """
	
	prefix = parameters[ parameters.index( '--prefix' )+1 ]
	if prefix[-1] != '/':
		prefix += '/'
	
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	seq_file = parameters[ parameters.index( '--cds' )+1 ]
	
	if not os.path.isfile( seq_file ):
		sys.exit( "ERROR: input CDS file not detected!" )
	
	if '--cpu' in parameters:
		cpu = int( parameters[ parameters.index( '--cpu' )+1 ] )
	else:
		cpu = 8
	
	RBH_file = prefix + "RBH_file.txt"
	blastdb = prefix + "blastdb"
	blast_result_file = prefix + "blast_result_file.txt"
	
	# --- running BLAST --- #
	if not os.path.isfile( blast_result_file ):
		os.popen( "/vol/biotools/bin/makeblastdb -in " + seq_file + " -out " + blastdb + " -dbtype nucl -parse_seqids" )
		os.popen( "/vol/biotools/bin/blastn -query " + seq_file + " -db " + blastdb + " -out " + blast_result_file + " -outfmt 6 -evalue 0.0001 -num_threads " + str( cpu ) )
	
	# --- loading data --- #
	seqs = load_multiple_fasta_file( seq_file )
	print "number of input sequences: " + str( len( seqs.keys() ) )
	data = load_results_from_BLAST_result_file( blast_result_file )
	RBHs = identify_RBHs( data, RBH_file )


if '--prefix' in sys.argv and '--cds' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
