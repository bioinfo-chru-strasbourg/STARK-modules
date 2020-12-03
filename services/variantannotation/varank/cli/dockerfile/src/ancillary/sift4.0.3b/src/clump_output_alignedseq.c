/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */ 

/*   clump_output_aligned.c
same as clump.c except output as aligned sequences
didn't add option to clump.c because would lose compatibility with other programs

clumps sequences (% clustering) and gives a consensus sequence for each cluster

most of this code is taken from Jorja's blweight.c

option 0 -- print out all of sequences (including those with only one sequence
	in cluster)
option 1 -- print out sequences that have 2 or more sequences in cluster

*/

#define EXTERN
#include "blocksprogs.h"
#include "PN_convert.c"
#include "Alignment.c"
#include "Matrix_Info.c"
#include "Psiblast.c"
#include "Clumping.c"

#define MINSEQ 2 
#define MAXSEQ 400
#define INDEX(n, col, row) (col*n - (col*(col+3))/2 - 1 + row)
#define MIN_CLUSTERS 3

void getargs (int argc, char* argv[], FILE** seqfp, 
		char outfilename[LARGE_BUFF_LENGTH], double* clus,
		int * option);

FILE* errorfp;
char errorfilename[LARGE_BUFF_LENGTH];

/* MAIN */

int 
main (int argc, char* argv[])
{
	FILE* seqfp; FILE* outfp; FILE* queryseqfp; FILE* keyoutfp;
	double clus;
	int nseqs, db_type, seq_type;
	Sequence* seqs[MAXSEQ];
	int no_of_clusters;
	int option;
	char outfilename[LARGE_BUFF_LENGTH];
	Sequence* queryseq;
	Sequence* consensus_seqs[MAXSEQ];
	int i;
	
	init_frq_qij();

        getargs (argc, argv, &seqfp, outfilename, &clus, &option);
        nseqs = 0;


     /*-----------------------------------------------------------------*/
      /*   Check next for input file of sequences & assume are aligned */
      db_type = type_dbs(seqfp, DbInfo);
      /*   could set db_type = FLAT if it comes back negative   */
      seq_type = UNKNOWN_SEQ;
      seq_type = seq_type_dbs(seqfp, DbInfo, db_type, seq_type);
      if (seq_type == NA_SEQ)
      {
         printf("WARNING: Sequences appear to be DNA but will be treated");
         printf(" as protein\n");
         seq_type = AA_SEQ;
      }
      rewind(seqfp);
      /*-----------------------------------------------------------------*/
      /*   read fasta sequences into memory                    */
      if (db_type >= 0)
      {
         while ( nseqs < MAXSEQ &&
             (seqs[nseqs] = read_a_sequence(seqfp, db_type, seq_type)) != NULL)
         {
            nseqs++;
         }
      }

       if ( (outfp = fopen (outfilename, "w")) == NULL)
        {
                printf ("Cannot open file %s \n", outfilename);
                exit (-1);
        }
        strcat (outfilename, ".consensuskey");
        if ( (keyoutfp = fopen (outfilename, "w")) == NULL)
        {
                printf ("Cannot open file %s\n", outfilename);
                exit (-1);
        }


	if (nseqs < MINSEQ) {
/*		printf ("There are %d sequences to find the motif in\n");
		printf ("Clumping automatically assigned to 1.0\n"); */
		clus = 1.0;
	} 
	
	/* remove sequences that have less than clus% 
			of length (should remove fragments)*/
	no_of_clusters = print_consensus_seqs_aligned (clus, seqs, nseqs, outfp, keyoutfp, option); 

	if (no_of_clusters < MIN_CLUSTERS) {
		fprintf (errorfp, "<BR>ERROR: After clustering at %d%%, only %d cluster(s) remained.  This suggests that there is not enough diversity among", 
					(int) (clus * 100), no_of_clusters);
		fprintf (errorfp, " the sequences.<BR>");
	       fprintf(errorfp, "*** We suggest you run a PSI-BLAST search on your ");
                fprintf (errorfp, "own and hand-pick sequences that");
                fprintf (errorfp, "may be related to your sequence.<BR>");
                fprintf (errorfp, "***You can then <A HREF=\"http://blocks.fhcrc.org/~pauline/SIFT_related_seqs_submit.html\">submit your ");
                fprintf (errorfp, "query sequence with the related sequences</A>");
                fprintf (errorfp, " or <BR><A HREF=\"http://blocks.fhcrc.org/~pauline/SIFT_aligned_seqs_submit.html\">submit the alignment.</A><BR>\n");
                fclose (errorfp);
                exit (-1);
	}

	fclose (errorfp);
	fclose (seqfp);
	fclose (outfp); fclose (keyoutfp);
	free_seqs (seqs, nseqs);
	rm_file (errorfilename);

	exit (0);

} /* end MAIN */

void getargs (int argc, char* argv[], FILE** seqfp, 
		char outfilename[LARGE_BUFF_LENGTH], double* clus,
		int* option)
{
        char seqfilename[LARGE_BUFF_LENGTH];
	char queryfilename[LARGE_BUFF_LENGTH];

	if (argc < 4) {
		printf ("clump.c : Clusters sequences into % clus");
	}

	if (argc > 1 ) strcpy (seqfilename, argv[1]);
        else
        {
                printf ("Enter filename with sequences:\n");
                gets (seqfilename);
        }
        if ((*seqfp = fopen (seqfilename, "r")) == NULL)
        {
                printf ("cannot open file %s \n", seqfilename);
                exit (-1);
        }

        if (argc > 2) strcpy (outfilename, argv[2]);
        else
        {
                printf ("Enter name of outfile\n");
                gets (outfilename);
        }


        strcpy (errorfilename, outfilename);
        strcat (errorfilename, ".error");
        if ((errorfp = fopen (errorfilename, "w")) == NULL) {
                printf ("couldn't open file %s\n", errorfilename);
                exit (-1);
        }


	if (argc > 3) *clus = atof (argv[3]);
	else
	{
		printf ("Enter clustering fraction (0.0 - 1.0)\n");
		scanf ("%lf", clus);
	}

	if (argc > 4) *option = atoi (argv[4]);
	else 
	{
		printf ("Enter option whether to print out clumps with just 1 seq\n");
		scanf ("%d", option);
	}

}

