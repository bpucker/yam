### Boas Pucker ###
### v0.3 ###
### bpucker@cebitec.uni-bielefeld.de ###

#updated with Araport11 annotation (March 2017)

__usage__ = """"
			python map_annotation_via_orthologs.py\n
			--in <INPUT_FILE>
			--out <OUTPUT_FILE>
			
			bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
			"""

import re, sys

# --- end of imports --- #

def load_annotation( annotation_file ):
	"""! @brief load all content from annotation file """
	
	mapping_table = {}
	with open( annotation_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if len( parts ) > 1:
				if parts[2] == "":
					annotation = "N/A"
				else:
					annotation = parts[1] + "." + parts[2]
				mapping_table.update( { parts[0].upper(): annotation } )
			line = f.readline()
	return mapping_table


def annotate_genes( mapping_table, data_input_file, output_file ):
	"""! @brief reads data from input file and maps annotation to geneID; results are written into output file """
	
	counter = 0
	with open( output_file, "w" ) as out:
		with open( data_input_file, "r" ) as f:
			f.readline()	#header
			line = f.readline()
			while line:
				parts = line.strip().split('\t')
				AGIs = re.findall( "AT[12345CM]{1}G\d{5}", line.upper() )
				if ',' in parts[1]:
					genes = parts[1].split(', ')
				else:
					genes = [ parts[1] ]
				anno = []
				for AGI in AGIs:
					try:
						anno.append( mapping_table[ AGI ] )
						counter += len( genes )
					except KeyError:
						anno.append( AGI + ". n/a" )
				anno = "_%%%_".join( anno )
				for gene in genes:
					new_line = [ gene, anno ]
					out.write( '\t'.join( new_line ) + '\n' )
				line = f.readline()
	print "number of annotated genes: " + str( counter )


def main( arguments ):
	""""! @brief runs all functions for mapping of annotation """
	
	annotation_file = "araport11_annotation.txt"
	mapping_table = load_annotation( annotation_file )
	
	data_input_file = arguments[ arguments.index( '--in' ) + 1 ]
	output_file = arguments[ arguments.index( '--out' ) + 1 ]
	
	annotate_genes( mapping_table, data_input_file, output_file )
	

if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
