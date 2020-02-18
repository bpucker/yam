### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python assess_read_len_distr.py
					--in <INPUT_FILE(FASTQ/FASTA)>
					--out <OUTPUT_FOLDER>
					
					optional:
					--type <INPUT_FILETYPE(fastq|fasta)[fastq]>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import os, sys, gzip
import matplotlib.pyplot as plt

# ---- end of imports --- #

def load_fastq_read_lens( input_file ):
	"""! @brief load all read lenths from given FASTQ file """
	
	read_lens = []
	if input_file[-2:] in [ "gz", "GZ" ] or input_file[-4:] in [ "gzip", "GZIP" ]:
		with gzip.open( input_file, "rb" ) as f:
			line = f.readline()	#header
			while line:
				seq = f.readline()
				read_lens.append( len( seq.strip() ) )
				f.readline()	#waste1
				f.readline()	#waste2
				line = f.readline()	#header
	else:
		with open( input_file, "r" ) as f:
			line = f.readline()	#header
			while line:
				seq = f.readline()
				read_lens.append( len( seq.strip() ) )
				f.readline()	#waste1
				f.readline()	#waste2
				line = f.readline()	#header
	return read_lens


def load_fasta_read_lens( input_file ):
	"""! @brief load read lengths from given FASTA file """
	
	read_lens = []
	if input_file[-2:] in [ "gz", "GZ" ] or input_file[-4:] in [ "gzip", "GZIP" ]:
		with gzip.open( input_file, "rb" ) as f:
			line = f.readline()	#header
			while line:
				seq = f.readline()
				read_lens.append( len( seq.strip() ) )
				line = f.readline()	#header
	else:
		with open( input_file, "r" ) as f:
			line = f.readline()	#header
			while line:
				seq = f.readline()
				read_lens.append( len( seq.strip() ) )
				line = f.readline()	#header
	return read_lens


def calculate_stats( read_lens ):
	"""! @brief calculate stats based on list of read lengths """
	
	total = sum( read_lens )
	sorted_read_lens = sorted( read_lens )[::-1]
	counter = 0
	i = 0
	while counter < ( 0.5*total ):
		counter += sorted_read_lens[ i ]
		i += 1
	n50 = sorted_read_lens[ i ]
	
	return total, n50


def generate_figure( read_lens, n50, fig_file ):
	"""! @brief generate histogram of read length distribution """
	
	cutoff = n50 * 3
	sorted_values = sorted( read_lens )
	data_to_plot = []
	for val in sorted_values:
		data_to_plot.append( min( [ val, cutoff ] ) )
	
	fig, ax = plt.subplots()
	
	ax.hist( data_to_plot, bins=100 )
	ax.set_xlabel( "read length [bp]" )
	ax.set_ylabel( "number of reads" )
	
	ax.set_xlim( 0, cutoff )
	
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	
	fig.savefig( fig_file )


def main( arguments ):
	"""! @brief run everything """
	
	input_file = arguments[ arguments.index( '--in' )+1 ]
	output_folder = arguments[ arguments.index( '--out' )+1 ]
	
	if '--type' in arguments:
		filetype = arguments[ arguments.index( '--type' )+1 ]
		if filetype not in [ "fasta", "fastq" ]:
			filetype = "fastq"
	else:
		filetype = "fastq"
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	fig_file = output_folder + "read_len_distr.pdf"
	
	print "recognized input file type: " + filetype
	
	# --- load data and generate data file --- #
	data_file = output_folder + "data.txt"
	if not os.path.isfile( data_file ):
		if filetype == "fastq":
			read_lens = load_fastq_read_lens( input_file )
		else:
			read_lens = load_fasta_read_lens( input_file )
		with open( data_file, "w" ) as out:
			out.write( "\n".join( map( str, read_lens ) ) + '\n' )
	
	# --- generate figure --- #
	with open( data_file, "r" ) as f:
		read_lens = map( int, f.read().strip().split('\n') )
	
	total, n50 = calculate_stats( read_lens )
	print "total amount of sequence: " + str( total )
	print "N50: " + str( n50 )
	generate_figure( read_lens, n50, fig_file )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
