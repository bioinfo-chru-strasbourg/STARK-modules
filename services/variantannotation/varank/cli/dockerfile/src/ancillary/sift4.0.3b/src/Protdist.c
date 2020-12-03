/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */   

#ifndef _PROTDIST_C_
#define _PROTDIST_C_

#include "Protdist.h"
#define FALSE 0
#define TRUE 1

void
protdist_format_out (Sequence* seqs[MAXSEQS], int nseqs, FILE* obfp)
{
        int i;

        fprintf (obfp, "%d %d\n", nseqs, seqs[0]->length);

        make_name_PHYLIP_width (seqs, nseqs);
        for (i = 0; i < nseqs; i++) {
                output_sequence_nonfasta (seqs[i], obfp);
        }

} /* end of protdist_format_out*/

 void
output_sequence_nonfasta (seq, osfp)
  Sequence* seq;
  FILE* osfp;
{
        int k;

        fprintf (osfp, "%s \n", seq->name);
        for (k=0;k < seq->length; k++) {
                fprintf (osfp, "%c", aa_btoa[seq->sequence[k]]);
               if ( ((k+1)%60 == 0) && (k != seq->length - 1) ) {
                        fprintf(osfp, "\n");
                }
        }
        fprintf (osfp, "\n");
} /* end of output_sequence_nonfasta */

/* make_name_PHYLIP_width : inserts whitespace after sequence name to
give a name length of PHYLIP_WIDTH characters */

void
make_name_PHYLIP_width (Sequence* seqs[MAXSEQS], int nseqs)
{
       int length, i, j;

        for (j=0; j < nseqs; j++) {
                length = strlen (seqs[j]->name);
                if (length < PHYLIP_WIDTH) {
                        for (i=length; i < PHYLIP_WIDTH; i++) {
                                seqs[j]->name[i] = ' ';
                        }
                	seqs[j]->name[PHYLIP_WIDTH] = '\0';
		} else if (length > PHYLIP_WIDTH) {
                        seqs[j]->name[PHYLIP_WIDTH] = '\0';
                }

        } /* end of for j */
} /* end of make_name_PHYLIP_width */

void
sort_scores_increasing_order (struct score_pair* pairwise_distances, int num_sequences)
{
/* use insertion sort because expect sequences to be in somewhat increasing distance */
	int i, j;

	double temp_score;
	char temp_name[15];

	for (i = 1; i < num_sequences; i++) {
		temp_score = pairwise_distances[i].score;
		strcpy (temp_name, pairwise_distances[i].name);
		for (j = i; j > 0 && pairwise_distances[j-1].score > temp_score; j--) {
			/* left score (j-1th entry) is bigger than jth (rightmost) score
				assign jth score j-1th's values */
			pairwise_distances[j].score = pairwise_distances[j-1].score;
			strcpy (pairwise_distances[j].name, pairwise_distances[j-1].name);
		}
		/* where j now belongs, j is smaller than i */
		pairwise_distances[j].score = temp_score;
		strcpy (pairwise_distances[j].name, temp_name);
	} /* end of for num_sequences */

} /* end sort_scores_increasing_order*/

int
seqs_with_threshold_diameter (struct score_pair* pairwise_distances, int num_sequences, double threshold_diameter)
{
	int i;
	double diameter;

	diameter = calculate_prodom_diameter (pairwise_distances, num_sequences);
	i = num_sequences - 1;
	while (i > 0 && diameter > threshold_diameter) {
		diameter = calculate_prodom_diameter (pairwise_distances, i);
		i--;
	}
	return (i + 1) ;
	
} /* end of seqs_with_threshold_diameter*/

double
calculate_prodom_diameter (struct score_pair* pairwise_distances, int num_sequences)
/* calculates root mean square from the first sequence */
{
	int i;
	double sum_of_squares;

	sum_of_squares = 0.0;

	for (i = 0; i < num_sequences; i++) {
		sum_of_squares += pow (pairwise_distances[i].score, 2);
	}
	return ( pow(sum_of_squares, .5));
	
} /* end of calculate_prodom_diameter */

int
seq_with_max_dist (struct score_pair* pairwise_distances, int num_sequences)
{
	double max_dist;
	int i;
	int index_of_max_seq;

	max_dist = 0;

	for (i = 0; i < num_sequences; i++) {
		if (pairwise_distances[i].score > max_dist) {
			index_of_max_seq = i;
			max_dist = pairwise_distances[i].score;
		}
	}
	return index_of_max_seq;
} /* end of seq_with_max_dist */


struct score_pair*
read_protdist (FILE* fp, int* num_seq)
{
	int i;
	struct score_pair* pairwise_distance;
	char line[LARGE_BUFF_LENGTH];
	char temp_line[LARGE_BUFF_LENGTH];
	char* strptr;
	int nseq;

	fscanf (fp, "%d\n", num_seq);

	nseq = *num_seq;

	pairwise_distance = (struct score_pair*) calloc (nseq, sizeof (struct score_pair));

/*	fgets (line, LARGE_BUFF_LENGTH, fp);
	strptr = strtok (line, " \t\n\r\0");
	strcpy (pairwise_distance[0].name, strptr);
	printf ("recording pairwise distances for %s\n", strptr);
*/
	fscanf (fp, "%s", pairwise_distance[0].name);
	printf ("1rst name %s \n", pairwise_distance[0].name);
	i = 0;
	while ( i < nseq) { 
		fscanf (fp, "%lf ", &pairwise_distance[i].score);
		/*pairwise_distance[i].score = atof (strptr); */
		i++;
	}
	/* now read in the rest of the names in order to correspond to scores */
	for (i = 1; i < nseq ; i++) {
		if (fgets (line, LARGE_BUFF_LENGTH, fp) == NULL) {
			printf ("Error, not enough sequences?\n");
			exit (-1);
		}
		strcpy (temp_line,line);
		strptr = strtok (temp_line, " "); /* temp_line has null after
						first word*/
		strcpy (pairwise_distance[i].name, strptr); 
		/* for lines longer than LARGE_BUFF_LENGTH, read in the rest
		   of the line */
		while (strrchr (line, '\n') == NULL) {
			fgets (line, LARGE_BUFF_LENGTH, fp);
		}
	}
	return pairwise_distance;
} /* end read_protdist */

void
free_array (struct score_pair* array)
{
	free (array);
	array = NULL;
}

			
#endif	
