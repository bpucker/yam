### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python reduce_to_tig.py
					--in <INPUT_FILE>
					--out <OUTPUT_FILE>
					
					WARNING: contig names need to be unique!
					(not working for splitted contigs e.g. after medaka polishing)
					"""

import re, sys

def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( fasta_file ) as f:
		header = re.findall( "tig\d+", f.readline() )[0]
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = re.findall( "tig\d+", line )[0]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def main( arguments ):
	"""! @brief run everything """
	
	input_file = arguments[ arguments.index('--in')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]

	seqs = load_sequences( input_file )
	with open( output_file, "w" ) as out:
		for key in seqs.keys():
			out.write( '>' + key + '\n' + seqs[ key ] + '\n' )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
