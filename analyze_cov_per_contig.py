### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python analyze_cov_per_contig.py
					--cov <COVERAGE_FILE>
					--out <OUTPUT_FOLDER>
					"""

import sys, os
import numpy as np
import matplotlib.pyplot as plt

# --- end of imports --- #

def analyze_coverage( cov_file,  cov_summary_file):
	"""! @brief analyze coverage per contig and summarize in one file """
	
	with open( cov_summary_file, "w" ) as out:
		with open( cov_file, "r" ) as f:
			line = f.readline()
			header = line.split('\t')[0]
			covs = []
			while line:
				parts = line.strip().split('\t')
				if parts[0] != header:
					out.write( header + '\t' + str( np.mean( covs ) ) + '\t' + str( np.median( covs ) ) + '\t' + str( len( covs ) )  + '\n' )
					covs = []
					header = parts[0]
				covs.append( float( parts[-1] ) )
				line = f.readline()
			out.write( header + '\t' + str( np.mean( covs ) ) + '\t' + str( np.median( covs ) ) + '\t' + str( len( covs ) ) + '\n' )


def load_coverage_values( cov_summary_file ):
	"""! @brief load average coverage value per chromosome """
	
	mean_cov_per_contig = {}
	median_cov_per_contig = {}
	contig_lengths = {}
	with open( cov_summary_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			mean_cov_per_contig.update( { parts[0]: float( parts[1] ) } )
			median_cov_per_contig.update( { parts[0]: float( parts[2] ) } )
			contig_lengths.update( { parts[0]: int( parts[3] ) } )
			line = f.readline()
	return mean_cov_per_contig, median_cov_per_contig, contig_lengths


def generate_figure( mean_cov_values, median_cov_values, figfile ):
	"""! @brief generate histogram of coverage distribution """
	
	fig, ax = plt.subplots()
	ax.hist( mean_cov_values, bins=5000, color="lime", alpha=0.5, label="mean" )
	ax.hist( median_cov_values, bins=5000, color="magenta", alpha=0.5, label="median" )
	ax.set_xlim( 0,200 )
	
	ax.legend()
	
	ax.plot( [ 15, 15 ], [ 0,20 ], color="black" )
	ax.plot( [ 50, 50 ], [ 0, 20 ], color="black" )
	
	ax.set_xlabel( "average coverage per contig" )
	ax.set_ylabel( "number of contigs" )
	
	plt.subplots_adjust( left=0.1, bottom=0.1, top=0.99, right=0.99 )
	
	fig.savefig( figfile, dpi=300 )


def calculate_genome_size( mean_cov_per_contig, contig_lengths ):
	"""! @brief calculate genome size based on coverage """
	
	genome_size = 0
	for contig in contig_lengths.keys():
		if mean_cov_per_contig[ contig ] < 30:
			genome_size += 2*contig_lengths[ contig ]
		elif mean_cov_per_contig[ contig ] < 150:
			genome_size += contig_lengths[ contig ]
		else:
			genome_size += contig_lengths[ contig ] * ( mean_cov_per_contig[ contig ] /  55 )
	
	print "estimated haploid genome size: " + str( genome_size / 2000000.0 ) + " Mbp"


def main( arguments ):
	"""! @brief run everything """

	cov_file = arguments[ arguments.index('--cov')+1 ]
	output_directory = arguments[ arguments.index('--out')+1 ]
	
	if output_directory[-1] != "/":
		output_directory += "/"
	
	if not os.path.exists( output_directory ):
		os.makedirs( output_directory )
	
	cov_summary_file = output_directory + "cov_summary_file.txt"
	figfile = output_directory + "cov_summary.pdf"

	analyze_coverage( cov_file,  cov_summary_file)
	mean_cov_per_contig, median_cov_per_contig, contig_lengths = load_coverage_values( cov_summary_file )
	generate_figure( mean_cov_per_contig.values(), median_cov_per_contig.values(), figfile )
	calculate_genome_size( mean_cov_per_contig, contig_lengths )


if '--cov' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
