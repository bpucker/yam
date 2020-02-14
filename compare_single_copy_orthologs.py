### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###


__usage__ = """
					python compare_single_copy_orthologs.py
					--ortho <ORTHOGROUP_FILE>
					--pep1 <PEPTIDE_FILE1>
					--pep2 <PEPTIDE_FILE2>
					--out <OUTPUT_FOLDER>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import os, sys, glob
import matplotlib.pyplot as plt

# ---- end of imports --- #


def get_single_copy_orthologs( orthofinder_result_file ):
	"""! @brief get all single copy orthologs """
	
	single_copy_orhthologs = {}
	with open( orthofinder_result_file, "r" ) as f:
		f.readline()	#header
		line = f.readline()
		while line:
			if not "," in line:
				parts = line.strip().split('\t')
				single_copy_orhthologs.update( { parts[1]: parts[2] } )
			line = f.readline()
	return single_copy_orhthologs


def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip().split( " " )[0]
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:].split( " " )[0]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def get_match_lines( alignment_file ):
	"""! @brief get all match lines """
	
	lines = []
	with open( alignment_file, "r" ) as f:
		f.readline()	#CLUSTAL header
		f.readline()	#empty line
		line = f.readline()
		while line:
			f.readline()	#empty line
			for i in range( 2 ):	#this function could be extended to allow inspection of more sequences
				line = f.readline()
			if len( line ) >= 16+60:
				lines.append( line[ 16:16+60 ] )
				#print line[ 16:16+60 ]
			else:
				lines.append( line[ 16:len( line )-1 ] )
				#print line[ 16:len( line )-1 ]
			line = f.readline()
	merged_match_line = "".join( lines )
	return merged_match_line.count("*"), merged_match_line.count(":"), merged_match_line.count("."), float( len( merged_match_line ) )


def main( arguments ):
	"""! @brief runs everything """
	
	orthofinder_result_file = arguments[ arguments.index( '--ortho' )+1 ]
	pep_file1 = arguments[ arguments.index( '--pep1' )+1 ]
	pep_file2 = arguments[ arguments.index( '--pep2' )+1 ]

	output_folder = arguments[ arguments.index( '--out' )+1 ]
	data_dir = output_folder + "data_dir/"

	if not os.path.exists( data_dir ):
		os.makedirs( data_dir )

	# --- find single copy orthologs --- #
	single_copy_orthologs = get_single_copy_orthologs( orthofinder_result_file )
	print "number of identified single copy orthologs: " + str( len( single_copy_orthologs.keys() ) )
	pep1 = load_sequences( pep_file1 )
	pep2 = load_sequences( pep_file2 )

	# --- combine both sequences in one file --- #
	for ortholog in single_copy_orthologs.keys():
		output_file = data_dir + ortholog + ".fasta"
		if not os.path.isfile( output_file ):
			with open( output_file, "w" ) as out:
				out.write( '>' + ortholog + '\n' + pep1[ ortholog ] + '\n>' + single_copy_orthologs[ ortholog ] + '\n' + pep2[ single_copy_orthologs[ ortholog ] ] + '\n' )

	# --- run mafft --- #
	fasta_files = glob.glob( data_dir + "*.fasta" )
	print "number of detected FASTA files: " + str( len( fasta_files ) )
	for fasta in fasta_files:
		mafft_output_file = fasta + ".clustal"
		if not os.path.isfile( mafft_output_file ):
			os.popen( "mafft --clustalout " + fasta + " > " + mafft_output_file )

	# --- analyze alignments --- #
	alignment_files = glob.glob( data_dir + "*.fasta.clustal" )
	print "number of detected alignment (CLUSTAL) files: " + str( len( alignment_files ) )
	identity = []
	similarity = []
	high_similarity = []
	combi1 = []
	combi2 = []
	for align in alignment_files:
		identical, similar, very_similar, alignment_length = get_match_lines( align )
		try:
			identity.append( identical / alignment_length )
			similarity.append( similar / alignment_length )
			high_similarity.append( very_similar / alignment_length )
			combi1.append( ( identical+very_similar ) / alignment_length )
			combi2.append( ( identical+very_similar + similar ) / alignment_length )
		except ZeroDivisionError:
			pass


	fig_file = output_folder + "alignment_similarity.pdf"

	fig, ax = plt.subplots(  )
	ax.hist( combi2, color= "magenta", label="identity+high similarity+similarity", bins=100 )
	ax.hist( combi1, color= "blue", label="identity+high similarity", bins=100 )
	ax.hist( identity, color= "lime", label="identity", bins=100 )

	ax.legend()
	ax.set_ylabel( "number of sequence pairs" )
	ax.set_xlabel( "proportion of residues in alignment" )

	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)

	fig.savefig( fig_file )



if '--ortho' in sys.argv and '--pep1' in sys.argv and '--pep2' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
