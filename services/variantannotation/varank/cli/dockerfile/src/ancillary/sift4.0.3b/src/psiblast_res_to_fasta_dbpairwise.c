/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */  

/* psiblast_res_to_fasta_dbpairwise.c

Take in psiblast results (pairwise output) and prints out results
of found sequences in fasta format.  

Output1:  All sequences found in fasta format.

Output2: Sequences found with residues for the X's in the first sequence
	 (query) removed so that the length of the alignment is the same
	 as the query sequence.

05-30-00 changed to read a filename of regions, to paste these regions 
	together
05-31-00 changed to read top-scoring, nonoverlapping pairwise alignment 
	in order to deal with multiple alignments for the same protein
06-15-00 allow to read when there's only one match

10-15-00 added error log for web.  errors automatically written to outfile.error

*/

#define EXTERN
#include <assert.h>
#include "blocksprogs.h"
#include "blocklist.c"
#include "Alignment.c"
#include "List_Number.c"
#include "Protdist.c"
#include "Psiblast.c"

#define MAXSEQ 400 /* Maximum number of sequences */
#define MIN_SEQ 5
#define LINE_LEN 800
#define FALSE 0
#define TRUE 1
#define MAXWIDTH 55 /* maximum block width */

FILE* errorfp; /* global variable, file is outfilename.error */
char errorfilename[LARGE_BUFF_LENGTH];

/* Local routines */
void getargs (int argc, char* argv[], FILE** psiblastfp, 
		 FILE** outfpqueryalign,
		 int* max_iterations, 
		char queryfilename[LARGE_BUFF_LENGTH]);

Sequence*  read_sequence_from_filename (char filename[LARGE_BUFF_LENGTH]);
/* ################	MAIN 	###########################  */

int main
(int argc, char* argv[])
{

	FILE* alignmentfp; FILE* regionsfp;
	FILE* outfpqueryalign; /* output file containing the alignment of the
				matches with the query sequence (gaps with
				query removed so that length of every sequence
				is the same as query) in fasta format */ 
	char queryfilename[LARGE_BUFF_LENGTH];
	int max_iterations;
	int converged;
	Sequence *fastaseqs[MAXSEQ];
	Sequence *alignedseqs[MAXSEQ];
	int nseqs, aa_length, i;
	Sequence* query_seq; Sequence* tmp_seq;
	int option_carve_gaps;	
	getargs (argc, argv, &alignmentfp, &outfpqueryalign,
		 &max_iterations, queryfilename);

/* READ QUERY SEQUENCE so that it can be included in the database to be made
   There is a limit to the # of FILE pointers, passing in the filename
   and doing it the stupid way by opening the file up in the function.*/
	query_seq = read_sequence_from_filename (queryfilename);

	query_seq->name[0] = '\0';
	query_seq->info[0] = '\0';
	strcpy (query_seq->name, "QUERY");

	printf ("query length %d\n", query_seq->length);
	nseqs = 0;

/*****READ PSIBLAST RESULTS ************************************/
	option_carve_gaps = TRUE; 
	/* sequences will be the amino acids that align to the query sequence*/
        
	read_psiblast_header_until_last (alignmentfp, max_iterations); 
	nseqs = psiblast_pairwise (alignmentfp, alignedseqs,
				   query_seq->length, option_carve_gaps); 
	fclose (alignmentfp); 
	alignmentfp = NULL;
	printf ("nseqs %d\n", nseqs);
	if (nseqs < MAXSEQ) { /* put query as first in sequence listed*/
			 	/* This is expected for a lot of the 
					subroutines */
		tmp_seq = alignedseqs[0];
		alignedseqs[0] = query_seq;
/*		copy_sequence (alignedseqs[0], query_seq); */
		alignedseqs[nseqs] = tmp_seq;
		nseqs++;  /* to count the query sequence */
	} 
	fix_names (nseqs, alignedseqs);
/******** WRITING OUT OUTPUT *************/
	for (i=0; i<nseqs; i++) {
       		convert_gap_to_X (alignedseqs[i]); 
		convert_sequence_with_beg_and_end_gaps_to_X (alignedseqs[i]);
		output_sequence (alignedseqs[i], outfpqueryalign);
	}
       	fclose (outfpqueryalign);

	free_seqs (alignedseqs, nseqs);
        if (nseqs < MIN_SEQ) {
                fprintf (errorfp, "ERROR! Not enough sequences (only %d) found by the PSI-BLAST search! Program terminated.\n", nseqs -1);
                exit (-1);
        } else if (nseqs > MAXSEQ) {
                fprintf (errorfp, "WARNING: exceeded maximum # of alignments in processing PSIBLAST result\n");
                exit (-1);
        }
	fclose (errorfp);
	rm_file (errorfilename);
	exit (0);


} /* end main */

void getargs (int argc, char* argv[], FILE** psiblastfp, 
		 FILE** outfpqueryalign,
                int* max_iterations, 
		char queryfilename[LARGE_BUFF_LENGTH]) 
{
	char psiblastfilename[LARGE_BUFF_LENGTH];
	char fasta_outfilename[LARGE_BUFF_LENGTH];
	char queryaligned_outfilename[LARGE_BUFF_LENGTH];
	
	if (argc < 3)
	{
		printf ("Psiblast results to fasta database \n");
		printf ("printing out matches in fasta format and");
		printf ("printing out the alignment to the query in");
		printf ("fasta format\n");
	}

	if (argc > 1) strcpy (psiblastfilename, argv[1]);
	else
	{
		printf ("Enter name of psiblast outfile with pairwise alignments:\n");
		gets (psiblastfilename);
	}

	if ((*psiblastfp = fopen (psiblastfilename, "r")) == NULL)
	{
		printf ("cannot open file %s \n", psiblastfilename);
		exit (-1);
	}

	if (argc > 2) strcpy (queryaligned_outfilename, argv[2]);
        else
        {
                printf ("Enter name of output file for which the psiblast \n");
		printf ("alignment to the query will be printed out in \n");
		printf ("fasta format\n");
                gets (queryaligned_outfilename);
        }

	if ((*outfpqueryalign = fopen (queryaligned_outfilename, "w")) == NULL)
	{
		printf ("cannot open file %s\n", queryaligned_outfilename);
		exit (-1);
	}
	strcpy (errorfilename, queryaligned_outfilename);
	strcat (errorfilename, ".error");
	if ((errorfp = fopen (errorfilename, "w")) == NULL) {
		printf ("couldn't open file %s\n", errorfilename);
		exit (-1);
	}

        if (argc > 3) *max_iterations = atoi (argv[3]);
        else {
                printf ("Enter the number of psiblast iterations\n");
                scanf ("%d", max_iterations);
        }

        if (argc > 4)  strcpy (queryfilename, argv[4]);
        else
        {
                printf ("Enter name of query file\n");
                scanf ("%s", queryfilename);
        }

 
} /* end of getargs */


Sequence*
read_sequence_from_filename (char filename[LARGE_BUFF_LENGTH])
{
	FILE* fp;
	Sequence* seq;

	if ((fp = fopen (filename, "r")) == NULL) {
		fprintf (errorfp, "couldn't open %s in processing PSIBLAST results\n", filename);
		exit (-1);
	}
	seq = read_a_sequence (fp, FASTA, AA_SEQ);
	fclose (fp);
	return seq;

} /* end of read_sequence_from_filename */
