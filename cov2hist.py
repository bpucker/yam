### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python cov2hist.py
					--cov <COVERAGE_FILE>
					--fig <FIGURE_FILENAME>
					"""

import sys, os
import numpy as np
import matplotlib.pyplot as plt

# --- end of imports --- #


def generate_figure( cov_values, figfile ):
	"""! @brief generate histogram of coverage distribution """
	
	fig, ax = plt.subplots()
	ax.hist( cov_values, bins=25000 )
	ax.set_xlim( 0,300 )
	
	ax.plot( [ 15, 15 ], [ 0,20 ], color="lime" )
	ax.plot( [ 50, 50 ], [ 0, 20 ], color="magenta" )
	
	ax.set_xlabel( "coverage" )
	ax.set_ylabel( "number of positions" )
	fig.savefig( figfile, dpi=300 )


def main( arguments ):
	
	cov_file = arguments[ arguments.index('--cov')+1 ]
	figfile = arguments[ arguments.index('--fig')+1 ]

	cov_values = []
	with open( cov_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			cov_values.append( int( float( parts[-1] ) ) )
			line = f.readline()

	generate_figure( cov_values, figfile )


if '--cov' in sys.argv and '--fig' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
