### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###

import sys, glob, re, os, time, datetime, shutil
import matplotlib.pyplot as plt
from operator import itemgetter

# --- end of imports --- #

__usage__ = """ python homology_based_gene_selection.py
							--query <INPUT_FILE (FASTA)>
							--subject <INPUT_FILE (FASTA)>
							--out <OUTPUT_FOLDER>
							feature requests and bug reports: bpucker@cebitec.uni-bielefeld.de
						"""


def submit_jobs_to_cluster( prefix, query_file_names, reference_blastp_db, para_jobs ):
	"""! @brief submit BLAST jobs for each file to cluster """
	
	IDs_to_check = []
	batch_ID = str( datetime.datetime.now() )[-3:]
	for idx, file_name in enumerate( query_file_names ):
		ID = "B_" + batch_ID + '_' + str( idx ).zfill(4)
		IDs_to_check.append( ID )
		sh_file = prefix + ID + '.sh'
		out_file = prefix + ID + '.out'
		err_file = prefix + ID + '.err'
		
		cmd = "blastp -query " + file_name + " -db " + reference_blastp_db + " -out " +  '.'.join( file_name.split('.')[:-1] ) + ".txt -outfmt 6 -evalue 0.00001"
		
		with open( sh_file, "w" ) as out:
				out.write( "#!/bin/bash\n" + " ".join( [ 	"echo " + '"',
																cmd + '"',
																"| qsub -cwd",
																"-N",
																ID,
																"-l vf=1G",
																"-l arch=lx-amd64",
																"-o",
																out_file,
																"-e",
																err_file
															] ) + '\n'
							)
		os.popen( "chmod +x " + sh_file )
		os.popen( sh_file )
		time.sleep(1)
		#os.remove( sh_file )
		waiting_status = True
		while waiting_status:
			qstat = os.popen( "qstat" )
			content = qstat.read()
			qstat_IDs = re.findall( "B_" + batch_ID + "_\d{4}", content )
			counter = 0
			for ID in qstat_IDs:
				if ID in IDs_to_check:
					counter += 1
			if counter < para_jobs:
				waiting_status = False
			else:
				time.sleep( 1 )
	
	waiting_status = True
	while waiting_status:
		qstat = os.popen( "qstat" )
		content = qstat.read()
		qstat_IDs = re.findall( "B_" + batch_ID + "_\d{4}", content )
		waiting_status = False
		for ID in IDs_to_check:
			if ID in qstat_IDs:
				for each in content.split('\n')[2:-1]:
					if ID in each.split()[2] and not 'd' in each.split()[4]:
						waiting_status = True
		time.sleep( 10 )


def load_sequences( multiple_fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( multiple_fasta_file ) as f:
		header = f.readline()[1:].strip()
		if " " in header:
			header = header.split(' ')[0]
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: seq } )
					header = line.strip()[1:]
					if " " in header:
						header = header.split(' ')[0]
					seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		sequences.update( { header: seq } )
	return sequences


def produce_multiple_query_files( query_file, cutoff, tmp_folder ):
	"""! @brief produce multiple query files """
	
	sequences = load_sequences( query_file )
	query_file_names = []
	
	len_counter = 0
	name_counter = 1
	query_file = tmp_folder + "0".zfill(4) + ".fasta"
	query_file_names.append( query_file )
	out = open( query_file, "w" )
	for idx, seq_id in enumerate( sorted( sequences.keys() ) ):
		if len_counter >= cutoff:
			len_counter = 0
			out.close()
			query_file = tmp_folder + str( name_counter ).zfill(4) + ".fasta"
			query_file_names.append( query_file )
			out = open( query_file, "w" )
			name_counter += 1
		out.write( '>' + seq_id + '\n' + sequences[ seq_id ] + '\n' )
		len_counter += len( sequences[ seq_id ] )
	out.close()
	return query_file_names


def final_processing( tmp_folder, final_result_file ):
	"""! @brief blt processing of BLAT results for identification of best hit """
	
	result_file_names = sorted( glob.glob( tmp_folder + "*.txt" ) )
	print "number of detected result files: " + str( len( result_file_names ) )
	os.chdir( tmp_folder )
	os.popen( "cat *.txt > " + final_result_file )
	print final_result_file


def load_best_blast_hit( blast_result_file ):
	"""! @brief load best blast hit per query """
	
	best_hits = {}
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				if float( parts[-1] ) > best_hits[ parts[0] ]['score']:
					best_hits[ parts[0] ] = { 'score': float( parts[-1] ), 'sim': float( parts[2] ), 'len': int( parts[3] ), 'subject': parts[1] }
			except:
				best_hits.update( { parts[0]: { 'score': float( parts[-1] ), 'sim': float( parts[2] ), 'len': int( parts[3] ), 'subject': parts[1] } } )
			line = f.readline()
	return best_hits


def rank_genes_by_score_ratio( self_hits, db_hits, output_file, seqs, score_ratio_figure ):
	"""! @brief rank genes based on score ratio between db and self """
	
	# --- calculate BLAST result statistics --- #
	data = []
	for key in sorted( seqs.keys() ):
		try:
			self = self_hits[ key ]
			try:
				other = db_hits[ key ]
				ratio = other['score'] / self['score']
				data.append( { 'ID': key, 'ratio': ratio, 'sim': other['sim'], 'len': other['sim'], 'coverage': float( other['len'] ) / len( seqs[ key ] ), 'self': "yes" } )
			except KeyError:
				data.append( { 'ID': key, 'ratio': 0, 'sim': 0, 'len': 0, 'coverage': 0, 'self': "yes" } )
		except KeyError:
			data.append( { 'ID': key, 'ratio': 0, 'sim': 0, 'len': 0, 'coverage': 0, 'self': "no" } )
	
	
	# --- write results into output file --- #
	data = sorted( data, key=itemgetter( 'ratio', 'coverage', 'sim', 'self' ) )[::-1]
	with open( output_file, "w" ) as out:
		for entry in data:
			out.write( "\t".join( map( str, [ entry['ID'], entry['ratio'], entry['coverage'], entry['sim'], entry['len'], entry['self'] ] ) ) + '\n' )
	
	
	# --- generate cumsum plot --- #
	x = []
	y = []
	for idx, each in enumerate( data ):
		x.append( idx )
		if idx == 0:
			y.append( each['ratio'] )
		else:
			y.append( y[-1] + each['ratio'] )
	
	fig, ax = plt.subplots()
	ax.plot( x, y, '.b-', color="lime" )
	ax.set_xlabel( "number of genes" )
	ax.set_ylabel( "cumulative score ratio" )
	fig.savefig( score_ratio_figure, dpi=600 )
	plt.close( "all" )


def main( arguments ):
	"""! @brief check inputs and call functions """
	
	if '--para_jobs' in arguments:
		para_jobs = int( arguments[ arguments.index( '--para_jobs' ) + 1 ] )
	else:
		para_jobs = 200
	
	query_file = arguments[ arguments.index( '--query' ) + 1 ]
	
	if '--splitting_cutoff' in arguments:
		cutoff = int( arguments[ arguments.index( '--splitting_cutoff' ) + 1 ] )
	else:
		cutoff=50000
	
	if '--cpus' in arguments:
		cpus = int( arguments[ arguments.index( '--cpus' ) + 1 ] )
	else:
		cpus = 20
	
	output_folder = arguments[ arguments.index( '--out' ) + 1 ]
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	if not '--db' in arguments:
		reference_blastp_db = output_folder + "ref_blast_db"
		subject_file = arguments[ arguments.index( '--subject' ) + 1 ]
		cmd = "makeblastdb -in " + subject_file + " -out " + reference_blastp_db + " -dbtype 'prot'"
		os.popen( cmd )
	else:
		reference_blastp_db = arguments[ arguments.index( '--db' ) + 1 ]
	
	# --- run BLASTp jobs on cluster against database --- #
	final_result_file = output_folder + "FINAL_RESULT_FILE_vs_DB.txt"
	if not os.path.isfile( final_result_file ):
		tmp_folder = output_folder + "tmp/"
		if not os.path.exists( tmp_folder ):
			os.makedirs( tmp_folder )
		query_file_names = produce_multiple_query_files( query_file, cutoff, tmp_folder )	#seq length increased compared to BLAT
		submit_jobs_to_cluster( tmp_folder, query_file_names, reference_blastp_db, para_jobs )
		final_processing( tmp_folder, final_result_file )
	
	self_final_result_file = output_folder + "FINAL_RESULT_FILE_vs_SELF.txt"
	if not os.path.isfile( self_final_result_file ):
		selfdb = output_folder + "selfdb"
		os.popen( "makeblastdb -in " + query_file + " -out " + selfdb + " -dbtype 'prot'" )
		os.popen( "blastp -query " + query_file + " -db " + selfdb + " -out " + self_final_result_file + " -outfmt 6 -evalue 0.00001 -num_threads " + str( cpus ) )
	
	db_hits = load_best_blast_hit( final_result_file )
	self_hits = load_best_blast_hit( self_final_result_file )
	
	score_ratio_file = output_folder + "SCORE_RATIO.txt"
	score_ratio_figure = output_folder + "SCORE_RATIO.pdf"
	
	seqs = load_sequences( query_file )
	rank_genes_by_score_ratio( self_hits, db_hits, score_ratio_file, seqs, score_ratio_figure )


if __name__ == '__main__':
	
	if '--query' in sys.argv and '--subject' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	elif '--query' in sys.argv and '--db' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
