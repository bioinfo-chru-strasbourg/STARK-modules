/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */ 
/* seqs_from_psiblast_res.c 

07-06-00 reads in sequences from psiblast results.  Takes sequences from a file
and prints out the complete sequence from psiblast for these sequences.

07-07-00 Takes the sequences and aligns them by psiblast flat master so
	no need to be processed by clustal

07-12-00 changed so it prints out entire sequence without gaps, needs to
	be processed by clustal
*/

#define EXTERN
#include <assert.h>
#include "blocksprogs.h"
#include "blocklist.c"
#include "Alignment.c"
#include "List_Number.c"
#include "Protdist.c"
#include "stringhash.c"
#include "Psiblast.c"

#define MAXSEQ 400 /* Maximum number of sequences */
#define LINE_LEN 800
#define FALSE 0
#define TRUE 1
#define MAXWIDTH 55 /* maximum block width */
#define MIN_SEQ 1

/* Local routines */
void getargs (int argc, char* argv[], FILE** psiblastfp, 
		 FILE** seqfp,
		 int* max_iterations, 
		char queryfilename[LARGE_BUFF_LENGTH],
		char outfilename[LARGE_BUFF_LENGTH], 
		char pid[SMALL_BUFF_LENGTH]);

Sequence*  read_sequence_from_filename (char filename[LARGE_BUFF_LENGTH]);

FILE* errorfp;
char errorfilename[LARGE_BUFF_LENGTH];

/* ################	MAIN 	###########################  */

int main
(int argc, char* argv[])
{

	FILE* psiblastfp; FILE* seqfp; FILE* outfp;
	FILE* outfpqueryalign; /* output file containing the alignment of the
				matches with the query sequence (gaps with
				query removed so that length of every sequence
				is the same as query) in fasta format */ 
	char queryfilename[LARGE_BUFF_LENGTH];
	char outfilename[LARGE_BUFF_LENGTH];
	char pid[SMALL_BUFF_LENGTH];
	int max_iterations; char* strptr;
	char line[LARGE_BUFF_LENGTH]; int converged;
	Sequence *seqs[MAXSEQ]; Sequence* newseqs[MAXSEQ];
	int nseqs, aa_length, i;
	Sequence* query_seq; Sequence* tmp_seq;
	HashTable seqnamehash;
	int option_carve_gaps;	
	char tmpfilename[LARGE_BUFF_LENGTH];
	char command_line[LARGE_BUFF_LENGTH];
	char temp[SMALL_BUFF_LENGTH];

	getargs (argc, argv, &psiblastfp, &seqfp,
		 &max_iterations, queryfilename, outfilename, pid);

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
        read_psiblast_header_until_last (psiblastfp, max_iterations); 
	nseqs = psiblast_pairwise (psiblastfp, seqs,
				   query_seq->length, option_carve_gaps); 
	fclose (psiblastfp); 
	
	if (nseqs > MAXSEQ) { 
		fprintf (errorfp, "WARNING: exceeded maximum # of alignments when getting sequences from psiblast results \n");
		fclose (errorfp);
		exit (-1);
	}
	fix_names (nseqs, seqs);

/*****STORE NAMES OF SEQUENCES OF INTEREST*********/
	seqnamehash = InitializeTable (30);
	while (!feof (seqfp) && fgets (line, LARGE_BUFF_LENGTH, seqfp) != NULL){
		strptr = strtok (line, " \t\n\0\r");
		strptr[strlen(strptr)] = '\0';
		Insert (strptr, seqnamehash);
		printf ("inserting %s\n", strptr);
	}
	fclose (seqfp);
	if (seqnamehash->no_of_elements < MIN_SEQ) {
		fprintf (errorfp, "%d sequence(s) were chosen. Less than the minimum number of sequences required.  \n", seqnamehash->no_of_elements);
		exit (-1);
	}

        if ((outfp = fopen (outfilename, "w")) == NULL)
        {
                fprintf (errorfp, "cannot open file %s when getting sequences from psiblastresults\n", outfilename);
                exit (-1);
        }
	output_sequence(query_seq, outfp);
	for (i=0; i<nseqs; i++) {
       	printf ("doe %s exist in sequname hash?\n", seqs[i]->name);
 	if (Exists (seqs[i]->name, seqnamehash) ) {
printf ("yes %s it does\n", seqs[i]->name);
			convert_sequence_with_beg_and_end_gaps_to_X(seqs[i]);
			output_sequence (seqs[i], outfp); 
		}
	}
       	fclose (outfp);

	free_seqs (seqs, nseqs); 
	DestroyTable(seqnamehash);
	rm_file (errorfilename);
	exit (0);

} /* end main */

void getargs (int argc, char* argv[], FILE** psiblastfp, 
		 FILE** seqfp,
                int* max_iterations, char queryfilename[LARGE_BUFF_LENGTH], 
		 char outfilename[LARGE_BUFF_LENGTH], char pid[SMALL_BUFF_LENGTH]) 
{
	char psiblastfilename[LARGE_BUFF_LENGTH];
	char seqfilename[LARGE_BUFF_LENGTH];
	
	if (argc < 5)
	{
		printf ("seqs_from_psiblast_res \n");
		printf ("reads in sequences from psiblast results.  ");
		printf ("Takes sequences from a file and prints out\n");
		printf ("the complete sequence from psiblast for these ");
		printf ("sequences.\n");
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

	if (argc > 2) strcpy (seqfilename, argv[2]);
        else
        {
                printf ("Enter name of seq file for which the psiblast \n");
		printf ("sequences will be printed \n");
                gets (seqfilename);
        }

	if ((*seqfp = fopen (seqfilename, "r")) == NULL)
	{
		printf ("cannot open file %s\n", seqfilename);
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

	if (argc > 5) strcpy (outfilename, argv[5]);
	else {
		printf ("Enter name of out file\n");
		scanf ("%s", outfilename);
	}

        strcpy (errorfilename, outfilename);
        strcat (errorfilename, ".error");
        if ((errorfp = fopen (errorfilename, "w")) == NULL) {
                printf ("couldn't open file %s\n", errorfilename);
                exit (-1);
        }


	if (argc > 6) strcpy (pid, argv[6]);
	else {
		printf ("Enter pid to be used as a suffix\n");
		scanf("%s", pid);
	} 
} /* end of getargs */


Sequence*
read_sequence_from_filename (char filename[LARGE_BUFF_LENGTH])
{
	FILE* fp;
	Sequence* seq;

	if ((fp = fopen (filename, "r")) == NULL) {
		printf ("couldn't open %s\n", filename);
		exit (-1);
	}
	seq = read_a_sequence (fp, FASTA, AA_SEQ);
	fclose (fp);
	return seq;

} /* end of read_sequence_from_filename */
