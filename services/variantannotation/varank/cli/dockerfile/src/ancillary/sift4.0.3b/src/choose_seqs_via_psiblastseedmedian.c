/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */           
/* choose_seqs_via_psiblast.c */
/* 02-14-01 use median threshold
05-18-01, changed to add PSIBLAST_TURN sequences at once, if threshold has not
	  been reached, just do another search.  If threshold has been reached,
	 recalculate median threshold for the 1 .. PSIBLAST_TURN sequences
	 that have been added.
	Should reduce time for calculation of median threshold O (length * seq)
	  by PSIBLAST_TURN x fold	 
*/
#define EXTERN
#include "blocksprogs.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "Matrix_Info.c"
#include "Protdist.c"
#include "stringhash.c"
#include "Psiblast.c" 
#include "PN_blocks.c"
#include "Clumping.c"

/* constants */
#define PSIBLAST_TURN 5

void getargs (int argc, char* argv[],
		 char query_seq_file[LARGE_BUFF_LENGTH],
		FILE** seqfp, 
		FILE** outfp,
		double* clumping_threshold, 
		char pid[SMALL_BUFF_LENGTH],
		double* median_threshold);

int write_seqs_not_in_hash (char query_database_file[LARGE_BUFF_LENGTH],
                         Sequence* seqs[MAXSEQ], int nseqs,
                                         HashTable hash);

void write_seq_to_file (char filename[LARGE_BUFF_LENGTH], Sequence* seq);

FILE* errorfp;
FILE* logfp;
char errorfilename[LARGE_BUFF_LENGTH];

/* MAIN */

int
main (int argc, char* argv[])
{
	FILE* outfp; FILE* temp_chk_fp;
	FILE* seqfp;
	int nseqs, db_type, seq_type;
	int i, k, index, nseqs_written; double dtemp;
	double clumping_threshold;
	double info_median;
	Block* workingblock; 
	HashTable seqnamehash;
	char seqname[KEY_WIDTH];
	char pid[SMALL_BUFF_LENGTH];
	char temp_chk_filename[LARGE_BUFF_LENGTH];
	char psiblastres_file[LARGE_BUFF_LENGTH];
	char query_database_file[LARGE_BUFF_LENGTH];
	char query_seq_file[LARGE_BUFF_LENGTH];
	char outfile[LARGE_BUFF_LENGTH];
	char string[10]; char query_filename[LARGE_BUFF_LENGTH];
	Sequence* seqs[MAXSEQ];
	Boolean done;
	double median;
	double median_threshold;
	struct ClusterOfSequences * cluster_of_seqs; int num_clusters;
	int num_seqs_removed;

	ErrorLevelReport = 5;

	init_frq_qij(); /* matrix conversion in order to calculate R */

	getargs (argc, argv, query_seq_file, &seqfp,  &outfp, 
		 &clumping_threshold, pid, &median_threshold);
	seqnamehash = InitializeTable (10);
	
	i = 0;
/**************** READ IN SEQUENCES (to get info comments) ***************/
	nseqs = 0;
	db_type = type_dbs (seqfp, DbInfo); seq_type = AA_SEQ;
	rewind (seqfp);
	while ( nseqs < MAXSEQ && 
		(seqs[nseqs] = read_a_sequence(seqfp, db_type, seq_type)) != NULL) {
		nseqs++;
	}
	fclose (seqfp);

printf ("clumping teshodl is %.2f\n",clumping_threshold);	
/**************** GET SEED BLOCK BY PERCENTAGE IDENTITY ******************/
/*       cluster_of_seqs = cluster_seqs (clumping_threshold, seqs, nseqs,
                                   &num_clusters);
	workingblock = make_block (seqs[0]->length, 0, cluster_of_seqs[0].nseqs, cluster_of_seqs[0].seqs, FALSE);

*/

/*only copy query sequence  (argumen 3 is 1 =>copy 1rst seq. only )*/
	workingblock = make_block (seqs[0]->length, 0, 1, seqs, FALSE); 
assert (workingblock != NULL);
	for (i = 0; i < workingblock->num_sequences; i++) {
		printf ("inserting %s in hash\n", workingblock->sequences[i].name);
		Insert (workingblock->sequences[i].name, seqnamehash);
	}
	median = calculate_median_information_of_block (workingblock, FALSE, TRUE); 

/* write the first sequence into a file that will be used as the query
for the psiblast search */
	strcpy (query_filename, pid);
        strcat (query_filename, ".TEMP_QUERY");
	write_seq_to_file (query_filename, seqs[0]);

assert (workingblock != NULL);
printf ("median is %.2f\n", median);
/* printf ("finished extracting\n"); */
/*	fprintf (outfp, "*********SEED BLOCK********\n"); 
	fprintf (outfp, "median %.3f \n", median);
	output_block_s (workingblock, outfp, FLOAT_OUTPUT);  
	fflush (outfp); */
	done = FALSE; 
/**********ITERATIVE LOOP *************************************/
while ( median > median_threshold && !done) { 
	/* before add next sequence, copy this block, which
	we know that passes criterion, to finalblock */


        /****PSIBLAST CALL  do it only every PSIBLAST_TURN round  */ 
		strcpy (temp_chk_filename, pid);
		strcat (temp_chk_filename, ".TEMP");
		strcpy (psiblastres_file, pid);
		strcat (psiblastres_file, ".TEMP_PSIBLAST");
		strcpy (query_database_file, pid);
		strcat (query_database_file, ".TEMP_DB");

/*strcpy (temp_chk_filename,"tmpout3pid.TEMP"); */
		nseqs_written = write_seqs_not_in_hash (query_database_file, 
					seqs, nseqs,
					 seqnamehash);
		if (nseqs_written == 0) {
		 	fprintf (stderr, "wrote 0 seqs to database\n");
			done = TRUE;
		} else {
		formatdb_system_call (query_database_file) ; 
		if ((temp_chk_fp = fopen (temp_chk_filename, "w")) == NULL) {
			fprintf (errorfp, "couldn't open checkpoint file %s when choosing sequences\n", temp_chk_filename);
			exit (-1);
		}
		block_to_checkpoint_file (workingblock, temp_chk_fp, logfp); 
		fclose (temp_chk_fp); 
     		printf ("prior to psiblast system call");
		psiblast_system_call (temp_chk_filename, query_database_file,
				      psiblastres_file, query_filename, logfp); 
		printf ("finished psiblast system call");
		} /* end of if seqs written to database file */
	/****END OF PSIBLAST CALL ****/

	/* GET TOP PSIBLAST_TURN HITS ******/
	for (k = 0; k < PSIBLAST_TURN && !done; k++) { 
		get_top_seq (seqname, seqnamehash, psiblastres_file); 
		printf ( "top sequence is %s\n", seqname); 
		if (seqname[0] != '\0' && seqname[0] != ' ') {
		/**** ADD TOP SEQ to BLOCK *******************/
			index = get_index_from_seqs (seqname, seqs, nseqs);
			assert (workingblock != NULL);
			/* block not weighted */
			add_sequence_to_block (workingblock,seqs[index], FALSE);
		        if (workingblock->num_sequences == nseqs) {
				/* added all the possible sequences to the block 
				have chosen all of the sequences, done*/
				done = TRUE;
			}
			Insert (seqname, seqnamehash);   /* update hash */        
		} else {
                        fprintf (stderr, "done\n");
                        done = TRUE;
                }
	} /* done adding PSIBLAST_TURN sequences ***/

	median = calculate_median_information_of_block (workingblock, 
								FALSE, FALSE);	
 	fprintf (logfp, "added %s %s median:%.3f median %d seq total\n",
                 	workingblock->sequences[workingblock->num_sequences-1].name,
			seqs[index]->info, median, workingblock->num_sequences); 
                 	fflush (logfp); fflush (stderr);
			fflush (errorfp); fflush (outfp); 

	/* prepare to run psiblast on next iteration */
       	fprintf (stderr, "rm search files\n");
	rm_file (temp_chk_filename); 
      	rm_file (psiblastres_file);
	strcat (query_database_file, ".*");
	rm_file (query_database_file);

} /* end of while loop */
	/* this is the working_block median */
	if (median > median_threshold) { /* didn't have enough sequences to
						get below threshold */

		fprintf (stderr, "WARNING: The sequences selected have a median of %.2f which is not as diverse as the median threshold %.2f that was requested.<BR>\n ",
			median, median_threshold);
/* commented out 07/16/00 so web won't quit. printing error to stderr rather than
error file */
/*		print_block_ids (workingblock, outfp, TRUE); */
		/* should free memeory, clean up later */
/*		fclose (outfp); fclose (errorfp); fclose (logfp);
		exit (0);*/ 
	/* exit 0 else will generate an error in shell that calls this program*/
	}
	/* now of the last PSIBLAST_TURN sequences added, find out which one gives
	   approaches the right median */	
	num_seqs_removed = 0;
	while (median < median_threshold) {
		workingblock->num_sequences--;
		median = calculate_median_information_of_block (workingblock, 
							FALSE, FALSE);
		num_seqs_removed++;
	}

	assert (num_seqs_removed <= PSIBLAST_TURN);

	fprintf (logfp, "********FINAL BLOCK ****** ");
	fprintf (logfp, "last block : median %.3f with %d sequences\n", median,
			workingblock->num_sequences); 
/*	output_block_s (workingblock, outfp, FLOAT_OUTPUT); */ 


/*	print_block_sequences (workingblock, outfp); */ 
	print_block_ids (workingblock, outfp, TRUE); 

	workingblock->num_sequences+= num_seqs_removed; 
				/* increase it back to deallocate
					memory properly */	
	free_seqs (seqs, nseqs);
	free_block (workingblock);
	DestroyTable (seqnamehash);
	fprintf (logfp, "SUCCESSFUL\n");
	fclose (outfp);	
	fclose (errorfp);
	fclose (logfp);	
	rm_file (errorfilename);
	exit (0);
} /* end of main */

/* return the number of sequences written */
int
write_seqs_not_in_hash (char query_database_file[LARGE_BUFF_LENGTH], 
			 Sequence* seqs[MAXSEQ], int nseqs,
                                         HashTable hash)
{

	FILE* outfp;
	int i;
	int nseqs_written;

	nseqs_written = 0;
	if ( (outfp = fopen (query_database_file, "w")) == NULL) {
		fprintf (errorfp, "couldn't write to %s \n", 
			query_database_file);
		exit (-1);
	}

	for (i = 0; i < nseqs; i++) {
		if (!Exists (seqs[i]->name, hash)) {
			nseqs_written++;
			output_sequence (seqs[i], outfp);
		}
	}

	fclose (outfp);
	return nseqs_written;
 
} /* end write_seqs_not_in_block */
		
void
getargs (int argc, char* argv[], 
	char query_seq_file[LARGE_BUFF_LENGTH],
	FILE** seqfp, 
	FILE** outfp, 
	double* clumping_threshold, 
	char pid[SMALL_BUFF_LENGTH],
	double* median_threshold)
{
	char seqfile[LARGE_BUFF_LENGTH];
	char outfile[LARGE_BUFF_LENGTH];
	char logfile[LARGE_BUFF_LENGTH];

	if (argc < 6)
	{
		printf ("CHOOSE SEQS VIA PSIBLAST SEED MEDIAN \n");
	}

        if (argc > 1) strcpy (query_seq_file, argv[1]);
        else {
                printf ("Enter query sequence file\n");
                gets (query_seq_file);
        }



        if ( argc > 2) strcpy (seqfile, argv[2]);
        else
        {
                printf ("Enter sequence file \n");
                gets (seqfile);
        }
        if ( (*seqfp = fopen (seqfile, "r")) == NULL )
        {
                printf ("Cannot open %s \n", seqfile);
                exit (-1);
        }

	if ( argc > 3) strcpy (outfile, argv[3]);
	else 
	{
		printf ("Enter outfile \n");
		gets (outfile);
	}

        strcpy (errorfilename, outfile);
        strcat (errorfilename, ".error");
        if ((errorfp = fopen (errorfilename, "w")) == NULL) {
                printf ("couldn't open file %s\n", errorfilename);
                exit (-1);
        }
	if ( (*outfp = fopen (outfile, "w")) == NULL )
        {
                printf ("Cannot open %s \n", outfile);
                exit (-1);
        }

	strcpy (logfile, outfile);
	strcat (logfile, ".log");
	if ((logfp = fopen (logfile, "w")) == NULL) {
		printf ("couldn't open log file %s\n", logfile);
		exit (-1);
	}

	if (argc > 4) *clumping_threshold = atof (argv[4]);
	else 
	{
		printf ("Enter clumping threshold\n");
		scanf ("%lf\n", clumping_threshold);
		printf ("the clumping threshold is %.2f\n", *clumping_threshold);
	}

	if (argc > 5) strcpy (pid, argv[5]);
	else {
		printf ("Enter process identification number\n");
		gets (pid);
		printf ("gets %s\n", pid);
	}
	
	if (argc > 6) *median_threshold = atof (argv[6]);
	else {
		printf ("Enter median threshold\n");
		scanf ("%lf\n", median_threshold);
	        printf ("the median threshold is %.2f\n", *median_threshold);
	}


} /* end of getargs */ 

void
write_seq_to_file (char filename[LARGE_BUFF_LENGTH], Sequence* seq) 
{
	FILE* fp;

	if ((fp = fopen (filename, "w")) == NULL) {
		fprintf (errorfp, "couldn't write to %s\n", filename);
		exit (-1);
	}
	output_sequence (seq, fp);
	fclose (fp);
}
