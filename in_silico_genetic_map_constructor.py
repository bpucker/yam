### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.15 ###


__usage__ = """
				python in_silico_genetic_map_constructor.py
				--assembly <ASSEMBLY_FASTA_FILE>
				--ref <REFERENCE_FASTA_FILE>
				--out <OUTPUT_FOLDER>
				
				optional:
				--cluster <flag activates submission of jobs to compute cluster>
				
				bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
				"""

import os, sys, glob, re, time, datetime, shutil
from operator import itemgetter

# --- end of imports --- #


def load_sequences_with_order( fasta_file, len_cutoff ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	order = []
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip().split( " " )[0]
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					order.append( header )
					header = line.strip()[1:].split( " " )[0]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )
		order.append( header )
	return sequences, order


def chunk_generator( sequence, kmer_len ):
	"""! @brief generate chunks of given sequence """
	chunks = []
	for i in range( len( sequence ) / kmer_len ):
		chunks.append( sequence[ i*kmer_len : (i+1) * kmer_len ] )
	return chunks


def generate_in_silico_markers( in_silico_marker_file, reference_file, len_cutoff, kmer_len ):
	"""! @brief generate in silico markers for in silico genetic map """
	
	seqs, seq_order = load_sequences_with_order( reference_file, len_cutoff )
	marker_counter = 0
	with open( in_silico_marker_file, "w" ) as out:
		for header in seq_order:
			seq = seqs[ header ]
			chunks = chunk_generator( seq, kmer_len )
			for idx, chunk in enumerate( chunks ):
				out.write( '>' + header + "_%_" + str( idx+1 ) + '\n' + chunk + '\n' )
			marker_counter += len( chunks )
	print "number of generated in silico markers: " + str( marker_counter )


def load_BLAST_results( blast_result_file, match_sim_cutoff, match_len_cutoff, unique_score_ratio ):
	"""! @brief load BLAST hits """
	
	# --- load raw hits --- #
	blast_hits = {}
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if int( parts[3] ) >= match_len_cutoff:
				if float( parts[2] ) >= match_sim_cutoff:
					try:
						blast_hits[ parts[0] ].append( { 'LG': parts[0].split('_%_')[0], 'cm': float( parts[0].split('_%_')[1] ), 'chr': parts[1], 'pos': ( int( parts[8] ) + int( parts[9] ) ) / 2.0, 'score': float( parts[-1] ) } )
					except KeyError:
						blast_hits.update( { parts[0]: [ { 'LG': parts[0].split('_%_')[0], 'cm': float( parts[0].split('_%_')[1] ), 'chr': parts[1], 'pos': ( int( parts[8] ) + int( parts[9] ) ) / 2.0, 'score': float( parts[-1] ) } ] } )
			line = f.readline()
	
	# --- screen and clean --- #
	final_hits = []
	for hits in blast_hits.values():
		if len( hits ) == 1:
			final_hits.append( hits[0] )
		else:
			sorted_hits = sorted( hits, key=itemgetter( 'score' ) )
			if sorted_hits[-2]['score'] / sorted_hits[-1]['score'] <= unique_score_ratio:
				final_hits.append( sorted_hits[-1] )
	return final_hits


def generate_genetic_map_file( allmaps_input_file, marker_info ):
	"""! @brief generate ALLMAPS input file """
	
	markers = sorted( marker_info, key=itemgetter('LG', 'cm') )
	with open( allmaps_input_file, "w" ) as out:
		out.write( "ScaffoldID\tScaffoldPosition\tLinkageGroup\tGeneticPosition\n" )
		for marker in markers:
			out.write( "\t".join( map( str, [ marker['chr'], marker['pos'], marker['LG'], marker['cm'] ] ) ) + '\n' )


def submit_jobs_to_cluster( prefix, query_file_names, reference_blastn_db, para_jobs, additional_options ):
	"""! @brief submit BLAST jobs for each file to cluster """
	
	IDs_to_check = []
	batch_ID = str( datetime.datetime.now() )[-3:]
	for idx, file_name in enumerate( query_file_names ):
		ID = "B_" + batch_ID + '_' + str( idx ).zfill(4)
		IDs_to_check.append( ID )
		sh_file = prefix + ID + '.sh'
		out_file = prefix + ID + '.out'
		err_file = prefix + ID + '.err'
		
		cmd = "/vol/biotools/bin/blastn -query " + file_name + " -db " + reference_blastn_db + " -out " +  '.'.join( file_name.split('.')[:-1] ) + ".txt " + additional_options
		
		with open( sh_file, "w" ) as out:
				out.write( "#!/bin/bash\n" + " ".join( [ 	"echo " + '"',
																cmd + '"',
																"| qsub -cwd",
																"-N",
																ID,
																"-l vf=1G",
																"-l arch=lx-amd64",
																"-P denbi",
																"-o",
																out_file,
																"-e",
																err_file
															] ) + '\n'
							)
		os.popen( "chmod +x " + sh_file )
		os.popen( sh_file )
		time.sleep(1)
		os.remove( sh_file )
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
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: seq } )
					header = line.strip()[1:]
					seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		sequences.update( { header: seq } )
	return sequences


def produce_multiple_query_files( query_file, cutoff ):
	"""! @brief produce multiple query files """
	
	prefix = query_file + str( datetime.datetime.now() )[-5:]  + "/"
	os.makedirs( prefix )
	
	sequences = load_sequences( query_file )
	
	query_file_names = []
	
	len_counter = 0
	name_counter = 1
	query_file = prefix + "0".zfill(4) + ".fasta"
	query_file_names.append( query_file )
	out = open( query_file, "w" )
	for idx, seq_id in enumerate( sorted( sequences.keys() ) ):
		if len_counter >= cutoff:
			len_counter = 0
			out.close()
			query_file = prefix + str( name_counter ).zfill(4) + ".fasta"
			query_file_names.append( query_file )
			out = open( query_file, "w" )
			name_counter += 1
		out.write( '>' + seq_id + '\n' + sequences[ seq_id ] + '\n' )
		len_counter += len( sequences[ seq_id ] )
	out.close()
	return prefix, query_file_names


def final_processing( query_file_names, final_result_file ):
	"""! @brief blt processing of BLAT results for identification of best hit """
	
	result_file_names = []
	for filename in query_file_names:
		result_file_names.append( '.fasta'.join( filename.split('.fasta')[:-1] ) + '.txt' )
	
	cmd1 = "cat " + " ".join( result_file_names ) + " > " + final_result_file
	os.popen( cmd1 )


def run_blastn_on_cluster( query_file, reference_blastn_db, output_file, para_jobs ):
	"""! @brief check inputs and call functions """
	
	additional_options = ' -outfmt 6 -evalue 0.01 '
	cutoff=5000000
	
	prefix, query_file_names = produce_multiple_query_files( query_file, cutoff )
	
	submit_jobs_to_cluster( prefix, query_file_names, reference_blastn_db, para_jobs, additional_options )
	
	final_processing( query_file_names, output_file )


def load_genetic_marker_infos( allmaps_input_file ):
	"""! @brief load infos from given input file """
	
	infos = []
	with open( allmaps_input_file, "r" ) as f:
		f.readline()
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			infos.append( { 'chr': parts[0], 'pos': int( float( parts[1] ) ), 'LG': parts[2], 'cm': parts[3] } )
			line = f.readline()
	return infos


def main( arguments ):
	"""! @brief run everything """
	
	assembly_file = arguments[ arguments.index( '--assembly' )+1 ]
	reference_file = arguments[ arguments.index( '--ref' )+1 ]
	output_dir = arguments[ arguments.index( '--out' )+1 ]
	
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	if '--cluster' in arguments:
		cluster_status = True
	else:
		cluster_status = False
	
	
	len_cutoff = 100000	#minimal length of reference sequences to be considered during genetic marker generation
	kmer_len = 1000	#lenght of in silico marker sequences
	cpus = 1	#number of cores to use for BLAST
	word_size = 12	#word size for BLAST of in silico markers
	match_sim_cutoff = 0.7	#minimal similarity required to keep BLAST hit
	match_len_cutoff = kmer_len*0.7	#minimal alignment length required to keep BLAST hit
	unique_score_ratio = 0.9	#max ratio between second best and best BLAST hit score to keep BLAST hit
	para_jobs = 500	#number of BLAST jobs to run in parallel
	
	# --- generating marker sequences --- #
	in_silico_marker_file = output_dir + "in_silico_markers.fasta"
	if not os.path.isfile( in_silico_marker_file ):
		generate_in_silico_markers( in_silico_marker_file, reference_file, len_cutoff, kmer_len )
	
	
	# --- map marker to assembly --- #
	blast_result_file = output_dir + "blast_results.txt"
	blastdb = output_dir + "blastdb"
	if not os.path.isfile( blast_result_file ):
		os.popen( "makeblastdb -in " + assembly_file + " -out " + blastdb + " -dbtype nucl" )
		if cluster_status:
			print "running BLASTn jobs on cluster ..."
			run_blastn_on_cluster( in_silico_marker_file, blastdb, blast_result_file, para_jobs )
		else:
			print "running BLAST locally... "
			os.popen( "blastn -query " + in_silico_marker_file + " -db " + blastdb + " -out " + blast_result_file + " -outfmt 6 -evalue 0.001 -word_size " + str( word_size ) +" -num_threads " + str( cpus ) )
	
	
	# --- generate genetic map file for ALLMAPS --- #	
	allmaps_input_file = output_dir + "in_silico_genetic_map.txt"
	if not os.path.isfile( allmaps_input_file ):
		marker_info = load_BLAST_results( blast_result_file, match_sim_cutoff, match_len_cutoff, unique_score_ratio )
		print "number of mapped genetic markers: " + str( len( marker_info ) )
		generate_genetic_map_file( allmaps_input_file, marker_info )
	
	
	# --- generate genetic map in BED format --- #
	allmaps_bed_file = output_dir + "genetic_map.bed"
	weight_file = output_dir + "weights.txt"
	if not os.path.isfile( allmaps_bed_file ):
		infos = load_genetic_marker_infos( allmaps_input_file )
		with open( allmaps_bed_file, "w" ) as out:
			for idx, x in enumerate( infos ):
				if idx % 20:	#### REDUCE THE NUMBER OF GENETIC MARKER BY ONLY TAKING A FEW OF THEM ###
					out.write( "\t".join( [ x['chr'], str( x['pos'] ), str( x['pos']+1 ), "genetic_map-" + x['LG']+":"+x['cm'], x['chr'] + ":" + str( x['pos']+1) ] ) + '\n' )
		with open( weight_file, "w" ) as out:
			out.write( "genetic_map\t1\n" )
	
	
	# --- run ALLMAPS --- #
	os.chdir( output_dir )
	os.popen( "PYTHONPATH=/vol/biotools/lib/ALLMAPS/lib64/python3.7/site-packages" )	#setting up environment
	cmd = [ 	"python3.7 -m jcvi.assembly.allmaps path",
					"-w " + output_dir + "weights.txt",
					allmaps_bed_file,
					assembly_file
				]
	os.popen( " ".join( cmd ) )
	os.popen( "python3.7 -m jcvi.assembly.allmaps estimategaps " + output_dir + "allmaps_input.bed" )


if '--assembly' in sys.argv and '--ref' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

