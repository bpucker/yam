### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python cov_hist_per_contig.py
					--cov <COVERAGE_FILE>
					--out <OUTPUT_FOLDER>
					"""

import matplotlib.pyplot as plt
import sys, os

# --- end of imports --- #


def construct_figure( fig_file_folder, header, cov_values ):
	"""! @brief construct coverage figure per contig """
	
	fig, ax = plt.subplots()
	ax.hist( cov_values, bins=10000 )
	ax.set_xlim( 0, 200 )
	ax.set_xlabel( "coverage" )
	ax.set_ylabel( "number of positions" )
	ax.set_title( header + "   len=" + str( len( cov_values ) ) )
	fig.savefig( fig_file_folder + header + ".pdf", dpi=300 )
	
	plt.close( "all" )


def main( arguments ):
	
	cov_file = arguments[ arguments.index( '--cov' )+1 ]
	fig_file_folder = arguments[ arguments.index( '--out' )+1 ]
	
	if fig_file_folder[-1] != '/':
		fig_file_folder += "/"
	
	if not os.path.exists( fig_file_folder ):
		os.makedirs( fig_file_folder )
	
	with open( cov_file, "r" ) as f:
		line = f.readline()
		cov_values = []
		header = line.split('\t')[0]
		while line:
			parts = line.strip().split('\t')
			if parts[0] != header:
				construct_figure( fig_file_folder, header, cov_values )
				header = parts[0]
				cov_values = []
			cov_values.append( int( float( parts[-1] ) ) )
			line = f.readline()


if '--cov' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
