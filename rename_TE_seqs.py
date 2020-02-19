### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python rename_TE_seqs.py
					--in <INPUT_FILE>
					--out <OUTPUT_FILE>
					"""

import os, sys, re

# --- end of imports --- #

def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( fasta_file ) as f:
		line = f.readline()[1:].strip()
		header = "family-" + re.findall( "family\-\d+", line )[0].split('-')[1].zfill( 5 ) + "\t(" + re.findall( "RepeatScout Family Size = \d+", line )[0]  + ")"
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					try:
						header = "family-" + re.findall( "family\-\d+", line )[0].split('-')[1].zfill( 5 ) + "\t(" + re.findall( "RepeatScout Family Size = \d+", line )[0]  + ")"
					except:
						try:
							header = "family-" + re.findall( "family\-\d+", line )[0].split('-')[1].zfill( 5 ) + "\t(" + re.findall( "Recon Family Size = \d+", line )[0]  + ")"
						except:
							print line
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def main( arguments ):
	"""! @brief run everything """
	
	input_file = arguments[ arguments.index( '--in' )+1 ]
	output_file = arguments[ arguments.index( '--out' )+1 ]
	
	seqs = load_sequences( input_file )
	print "number of TE consensi: " + str( len( seqs.keys() ) )
	with open( output_file, "w" ) as out:
		for key in seqs.keys():
			out.write( '>' + key + '\n' + seqs[ key ] + '\n' )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
