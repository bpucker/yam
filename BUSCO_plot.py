### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python BUSCO_plot.py
					--in <INPUT_FILE_WITH_BUSCO_SUMMARY_STRINGS>
					--fig <FIGURE_FILE (PDF)>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import re, sys, os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

# --- end of imports --- #

def load_data( busco_summary_file ):
	"""! @brief load all data from given file """
	
	results = []
	with open( busco_summary_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				result_string = re.findall( "C:\d+\.\d+%\[S:\d+\.\d+%,D:\d+\.\d+%\],F:\d+\.\d+%,M:[-]*\d+\.\d+%,n:\d+", parts[1] )[0]
				results.append( { 	'id': parts[0],
											'c': float( result_string.split('C:')[-1].split('%')[0] ),
											's': float( result_string.split('S:')[-1].split('%')[0] ),
											'd': float( result_string.split('D:')[-1].split('%')[0] ),
											'f': float( result_string.split('F:')[-1].split('%')[0] ),
											'm': float( result_string.split('M:')[-1].split('%')[0] ),
											'n': int( result_string.split('n:')[-1] )
										} )
			except IndexError:
				print line
			line = f.readline()
	return results


def generate_plot( results, fig_file ):
	"""! @brief generate stacked bar plots """
	
	fig, ax = plt.subplots()
	
	values = []
	labels = []
	for each in results:
		values.append( [ each['s'], each['d'], each['f'], each['m'] ] )
		labels.append( each['id'] )
	
	my_plots = []
	colors = [ "green", "#4daf4a", "cyan", "grey" ]
	for idx, each in enumerate( values ):
		tmp_plot = []
		for i in range( 4 ):
			if i > 0:
				tmp_plot.append( ax.bar( idx, each[i], width=0.8, bottom=sum( each[ :i ] ), color=colors[i] ) )
				if each[i] < 4:
					ax.text( idx, sum( each[ :i ] )+each[i]*0.5, str( each[i] )+"%", ha="center", va="center", fontsize=6 )
				else:
					ax.text( idx, sum( each[ :i ] )+each[i]*0.5, str( each[i] )+"%", ha="center", va="center", fontsize=10 )
			else:
				tmp_plot.append( ax.bar( idx, each[i], width=0.8, color=colors[i] ) )
				ax.text( idx, each[i]*0.5, str( each[i] )+"%", ha="center", va="center" )
		my_plots.append( tmp_plot )


	my_legend = [ 	mpatches.Patch( color='grey', label='missing' ),
							mpatches.Patch( color='cyan', label='fragmented' ),
							mpatches.Patch( color='#4daf4a', label='duplicated' ),
							mpatches.Patch( color='green', label='single' )
						]
	ax.legend( handles=my_legend, fontsize=8, bbox_to_anchor=(0.96, 0.9) )	#loc=2, 


	#ax.set_xlabel( "assemblies", fontsize=12 )
	ax.set_ylabel( "percentage of detected BUSCOs", fontsize=10 )

	ax.spines['right'].set_color('none')
	ax.spines['top'].set_color('none')

	start, end = ax.get_xlim()
	ax.xaxis.set_ticks( np.arange( 0, len( labels ), 1) )

	ax.set_xticklabels( labels, fontsize=8 )

	plt.subplots_adjust( bottom=0.075, left=0.1, right=0.85, top=0.99 )
	
	fig.savefig( fig_file, dpi=900 )


def main( arguments ):
	"""! @brief run everything """
	
	busco_summary_file = arguments[ arguments.index('--in')+1 ]
	fig_file = arguments[ arguments.index('--fig')+1 ]
	
	results = load_data( busco_summary_file )
	generate_plot( results, fig_file )


if '--in' in sys.argv and '--fig' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
