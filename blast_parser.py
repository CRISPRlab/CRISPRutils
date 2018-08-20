#!/usr/bin/env python

# USAGE: to use on xml results from a cmd line BLAST using -outfmt 5
#		blast_parser.py [-f blast.xml - required] 
#						[-c/--clean (True/False) - optional] 
#						[-t/--trim (True/False) - optional (removes lines not matching inclusionText: used with -c)]
#						[-p/--predict (True/False) - optional] 
#						[-x/--remote results - to use when BLAST results were obtained against a remote database]
#						


from Bio.Blast import NCBIXML
import sys
from os import walk
from pprint import pprint
from Bio import Entrez
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
import argparse
import pandas
import os
import os.path
import Bio
import warnings

# parse input arguments #
parser = argparse.ArgumentParser(description='BlastParser by Matt Nethery')
parser.add_argument('-f','--inputFile', help='Input xml filename',required=True)
parser.add_argument('-c','--clean', dest='clean', action='store_true')
parser.set_defaults(clean=False)
parser.add_argument('-t','--trim', dest='trim', action='store_true')
parser.set_defaults(trim=False)
parser.add_argument('-p','--predict', dest='predict', action='store_true')
parser.set_defaults(predict=False)
parser.add_argument('-r','--rev_comp', dest='rev_comp', action='store_true')
parser.set_defaults(rev_comp=False)
parser.add_argument('-x','--remote', dest='remote', action='store_true')
parser.set_defaults(remote=False)
args = parser.parse_args()

# in some cases packages were compiled against an older numpy than is installed #
# mute these messages #
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
 
## show values ##
print ("Input filename: %s" % args.inputFile)
print ("Trim?: %s" % args.trim)
print ("Predict PAM?: %s" % args.predict)

#Modify inclusion text as needed
inclusionText = [
					"plasmid"
					, "phage"
					, "bacteriophage"
					, "prophage"
				]

## global vars ##
margin = 10


def errorAndExit() :
	print ("ERROR: You must provide an appropriate BLAST xml file result as the primary argument")
	exit()

if args.inputFile and len(args.inputFile) > 0 :
	result_handle = open(args.inputFile)
else :
	errorAndExit()

if args.clean :

	blast_records = NCBIXML.parse(result_handle)
	
	masterList = []
	#output tab delmited file
	if os.path.isfile(args.inputFile + '.tab.hits.csv') :
		os.remove(args.inputFile + '.tab.hits.csv')
	
	outputFile = open(args.inputFile + '.tab.hits.csv', 'w')
	outputFile.write("Query_Seq"+"\t"+"Accession"+"\t"+"Accession_Desc"+"\t"+"Expect"+"\t"+"Score"+"\t"+"Identity_Score"+"\t"+"Identity"+"\t"+"Positives"+"\t"+"Gaps"+"\t"+"Length"+"\t"+"Query_Start"+"\t"+"Match_Start"+"\t"+"Original_Query_Length"+"\t"+"Original_Identity_Score"+"\t"+"Strand"+"\n")
	for blast_record in blast_records:
	
		query = blast_record.query
		original_query_length = blast_record.query_length
	
		#hits:
		hits = []
		for alignment in blast_record.alignments :
			# tmpDict = {'match' : alignment.title}


			if args.remote:
				### for traditional remote blasts against nt/nr ###
				splitArr = alignment.title.split("|")
				accession = splitArr[3]
				hitDesc = splitArr[4]
				tmpDict = {'accession' : accession, 
							'accessionDesc' : hitDesc}
			else:
				### for local blasts against db's you built locally: ###
				tmpDict = {'accession' : alignment.title.split("|")[3], 
						'accessionDesc' : alignment.title}
					
	
			if args.trim :
				if not any(x in tmpDict['accessionDesc'] for x in inclusionText) :
					continue
	
			for hsp in alignment.hsps :
				tmpDict['score'] = hsp.score
				tmpDict['expect'] = hsp.expect
				tmpDict['gaps'] = hsp.gaps
				tmpDict['spcerSeq'] = hsp.sbjct
				tmpDict['matchSeq'] = hsp.query
				tmpDict['bits'] = hsp.bits
				tmpDict['identity'] = hsp.identities
				tmpDict['identityScore'] = str(round((float(hsp.identities) / float(hsp.align_length)) * 100, 2))
				tmpDict['length'] = hsp.align_length
				tmpDict['positives'] = hsp.positives
				tmpDict['frame'] = hsp.frame
				tmpDict['queryStart'] = hsp.sbjct_start #match
				tmpDict['subjectStart'] = hsp.query_start #spacer
				tmpDict['frame'] = hsp.frame
				tmpDict['originalQueryLength'] = original_query_length
				tmpDict['originalQueryScore'] = str(round((float(hsp.identities) / float(original_query_length)) * 100, 2))
	
				#add lines to output file
				outputFile.write(query+"\t"+tmpDict['accession']+"\t"+tmpDict['accessionDesc']+"\t"+str(float(tmpDict['expect']))+"\t"+str(tmpDict['score'])+"\t"+tmpDict['identityScore']+"\t"+str(tmpDict['identity'])+"\t"+str(tmpDict['positives'])+"\t"+str(tmpDict['gaps'])+"\t"+str(tmpDict['length'])+"\t"+str(tmpDict['subjectStart'])+"\t"+str(tmpDict['queryStart'])+"\t"+str(tmpDict['originalQueryLength'])+"\t"+str(tmpDict['originalQueryScore'])+"\t"+str(tmpDict['frame'])+"\n")
	
	
			if len(masterList) > 0 and masterList[-1]['query'] == query :
				lastObj = masterList[-1]
				hitList = lastObj['hits']
				hitList.append(tmpDict)
			else :
				hits.append(tmpDict)
				masterDict = {
								'query' : query, 
								'hits' : hits
								}
				masterList.append(masterDict)
		
	
	# pprint(masterList)
	outputFile.close()


if args.predict :

	#begin entrez seq downloads by accession #
	Entrez.email = ""
	Entrez.tool = "BlastParser.py"
	
	outputFile = open('pam_predict_spacer_list.fa', 'w')
	flankFile = open('pam_predict_flanks.fa', 'w')
	
	#import csv results and serialze
	fileRef = args.inputFile + '.tab.hits.csv'

	masterCsvData = pandas.read_csv(fileRef, index_col=None, sep='\t')
	masterCsvDict = [item for item in masterCsvData.T.to_dict().values()]
	queryList = [] # used to ensure duplicate spacers aren't searched (it misrepresents the PAM)
	for dict in masterCsvDict :
		query = dict['Query_Seq']
		print(query)
		if query in queryList: # prevent duplicates
			continue

		queryList.append(query)

		print('Accession: >>', dict['Accession'])
		tmpAccession = dict['Accession']

		handle = Entrez.efetch(db="nucleotide", id=tmpAccession, rettype="gb", retmode="text")
		record = SeqIO.read(handle, "genbank")
		handle.close()
		protoSeq = (record.seq)	

		length = dict['Length']
		expected_spacer_length = dict['Original_Query_Length']
		strand = str(dict['Strand']).strip("()")
		strandDir = int(strand.split(',')[1])	

		adjust = 0
		if length < expected_spacer_length :
			adjust = dict['Query_Start'] - 1	

		adjust = adjust if strandDir == 1 else -adjust

		protoStart = dict['Match_Start'] - adjust
		protoEnd = (protoStart + expected_spacer_length - 1) if strandDir == 1 else (protoStart - expected_spacer_length + 1)	
		if strandDir == 1 :
			protospacer = protoSeq[(protoStart-1-margin):protoEnd+margin]
		else :
			tmpSeq = protoSeq[protoEnd-1-margin:protoStart+margin]
			protospacer = Seq(str(tmpSeq), generic_dna).reverse_complement()			
		print (protospacer)

		outputFile.write(">" + query + "_" + dict['Accession'] + " >>" + dict['Accession_Desc'] + "\n")
		if args.rev_comp :
			protospacer = Seq(str(protospacer), generic_dna).reverse_complement()
		outputFile.write(str(protospacer) + "\n")

		upFlank = protospacer[0:10]
		downFlank = protospacer[-10:]
		flankFile.write(">" + query + "_" + dict['Accession'] + " >>" + dict['Accession_Desc'] + "\n")
		flankFile.write(str(upFlank+"--"+downFlank) + "\n")
		print (upFlank + " " + downFlank)

	outputFile.close()







