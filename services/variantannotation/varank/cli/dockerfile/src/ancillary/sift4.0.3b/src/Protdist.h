/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */   

#ifndef _PROTDIST_H_
#define _PROTDIST_H_

#define PHYLIP_WIDTH 10
#include "stringhash.c"

/* manipulates and reads protdist file.  
   routines for using protdist information 
*/

struct score_pair {
        char name[PHYLIP_WIDTH + 1];
        double score; /* record score query has with name */
};

void protdist_format_out (Sequence* seqs[MAXSEQS], int nseqs, FILE* obfp);
/* outputs sequences in sequential format for ProtDist */

void make_name_PHYLIP_width (Sequence* seqs[MAXSEQS], int nseqs);

void output_sequence_nonfasta ();
/* Sequence *seq, FILE* osfp .
   can pass in Sequence structure, and changes will still be made */

double calculate_prodom_diameter (struct score_pair* pairwise_distances, int num_sequences );

int seq_with_max_dist (struct score_pair* pairwise_distances, int num_sequences);

struct score_pair* read_protdist (FILE* fp, int* num_seq);

void sort_scores_increasing_order (struct score_pair* pairwise_distances, int num_sequences);

void free_array (struct score_pair* array);

int seqs_with_threshold_diameter (struct score_pair* pairwise_distances, int num_sequences, double threshold_diameter);


#endif /* protdist.h*/

