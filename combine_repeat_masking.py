### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python combine_repeat_masking.py
					--in1 <INPUT_FILE2>
					--in2 <INPUT_FILE1>
					--out <OUTPUT_FILE>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import sys

# --- end of imports --- #

def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ).upper() } )
					header = line.strip()[1:]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ).upper() } )	
	return sequences


def construct_new_seqs( seqs1, seqs2 ):
	
	new_seqs = {}
	for key in seqs1.keys():
		novel = []
		for idx, nt in enumerate( seqs1[ key ] ):
			if nt == "N":
				novel.append( "N" )
			elif seqs2[ key ][ idx ] == "N":
				novel.append( "N" )
			else:
				novel.append( nt )
		new_seqs.update( { key: "".join( novel ) } )
	return new_seqs


def main( arguments ):
	
	masked_file1 = arguments[ arguments.index( '--in1' )+1 ]
	masked_file2 = arguments[ arguments.index( '--in2' )+1 ]

	final_masked_file = arguments[ arguments.index( '--out' )+1 ]

	seqs1 = load_sequences( masked_file1 )
	seqs2 = load_sequences( masked_file2 )

	new_seqs = construct_new_seqs( seqs1, seqs2 )

	with open( final_masked_file, "w" ) as out:
		for key in new_seqs.keys():
			out.write( '>' + key + '\n' + new_seqs[ key ] + '\n' )


if '--in1' in sys.argv and '--in2' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
