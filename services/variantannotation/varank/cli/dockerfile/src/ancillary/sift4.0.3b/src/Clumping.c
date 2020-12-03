/* COPYRIGHT 2000 , Fred Hutchinson Cancer Research Center */ 

#ifndef _CLUMPING_C_
#define _CLUMPING_C_

/* 02-26-01 played around with cluster_seqs. didn't end up changing 
		anything, but
		didn't recompile to check 
*/
 
struct ClusterOfSequences {
	Sequence* seqs[MAXSEQ];
	int nseqs;
}; 

/* already declared in Alignment.c comment out 
struct cluster_pair {
        double score;
        int cluster;
};
*/

void free_ClusterOfSequences (struct ClusterOfSequences* cluster, 
				int num_clusters)
{
	int i;

/*	for (i = 0; i < num_clusters; i++) {
		free (cluster[i].seqs);
	} */
	free(cluster);
}


#define INDEX(n, col, row) (col*n - (col*(col+3))/2 - 1 + row)

struct ClusterOfSequences*
cluster_seqs (double clus, Sequence* seqs[MAXSEQ], int nseq, int* numclusters)
{

   int iclus, npair, threshold, s1, s2, l1, l2, px, i, i1, i2, first, j;
   int check1, check2, temp_cluster_no, tempint;
   int nclus[MAXSEQ], icluster[MAXSEQ], minclus, oldclus, width;
   int match, no_of_clumps;
   struct cluster_pair *pairs;
   int alignment_width;
  int nseqs_cluster, cluster_num;
  Sequence* seqs_cluster[MAXSEQ]; Block* block; Matrix* pssm;
        Sequence* consensus_sequence;
        char comment[SMALL_BUFF_LENGTH];
  struct ClusterOfSequences * clusters_of_seqs; 


   no_of_clumps = 0;
   npair = nseq*(nseq-1)/2;
   pairs = (struct cluster_pair *) malloc(npair * sizeof(struct cluster_pair));

   /*    Compute scores for all possible pairs of sequences            */
   for (s1=0; s1<nseq-1; s1++)                  /* col = 0, n-2     */
   {
      for (s2=s1+1; s2<nseq; s2++)      /* row = col+1, n-1 */
      {
         alignment_width = seqs[0]->length;
         match = 0;
         px = INDEX(nseq, s1, s2);
         pairs[px].score = 0;
         pairs[px].cluster = -1;
         for (i=0; i< seqs[0]->length; i++)
         {
                if ( seqs[s1]->sequence[i] == seqs[s2]->sequence[i]) {
                        if (seqs[s1]->sequence[i] == aa_atob['X'] ||
                              seqs[s1]->sequence[i]== aa_atob['-']) {
                                alignment_width--;
                        } else {
                                match++;
                        }
                } /* end of if letters match */
          } /* end of for width */
         pairs[px].score = ((double) match)/ ((double) alignment_width);
      }  /* end of s2 */
   }  /* end of s1 */

   /*-------Cluster if score exceeds threshold by scanning cols (s1) */
   for (s1=0; s1<nseq; s1++)
   {
      icluster[s1] = -1;                        /* clear out old values */
      nclus[s1] = 0;
   }
   iclus = 0;                                   /* cluster number */
   for (s1=0; s1<nseq-1; s1++)                  /* col = 0, n-2     */
      for (s2=s1+1; s2<nseq; s2++)      /* row = col+1, n-1 */
      {
         px = INDEX(nseq, s1, s2);
         if (pairs[px].score >= clus)      /*  cluster this pair */
         {
            if (icluster[s1] < 0)          /* s1 not yet clustered */
            {
               if (icluster[s2] < 0)       /* new cluster */
               {
                 icluster[s1] = iclus++;
                  icluster[s2] = icluster[s1];
                }
               else                             /* use s2's cluster  */
                  icluster[s1] =  icluster[s2];
            }
            /*  use s1's cluster if it has one and s2 doesn't */
            else if (icluster[s1] >= 0 && icluster[s2] < 0)
               icluster[s2] = icluster[s1];
            /* merge the two clusters into the lower number */
            else if (icluster[s1] >= 0 && icluster[s2] >= 0)
            {
               minclus = icluster[s1]; oldclus = icluster[s2];
               if (icluster[s2] < icluster[s1])
               {
                  minclus = icluster[s2]; oldclus = icluster[s1];
               }
               for (i1=0; i1<nseq; i1++)
                 if (icluster[i1] == oldclus)
                     icluster[i1] = minclus;
            }
         }  /* end of if pairs */
      }  /* end of s2 */

   /*---  Set ncluster, get rid of negative cluster numbers --*/
   for (s1=0; s1<nseq; s1++)
   {
      if (icluster[s1] < 0)
          icluster[s1] = iclus++;
   }
   for (s1=0; s1<nseq; s1++)
          nclus[icluster[s1]] += 1;

   i2 = 0;
 printf ("no of clusters %d\n", iclus); 


   /* move QUERY to first cluster */
        if (icluster[0] != 0) { /* QUERY is not first in cluster */
printf ("entered if statment\n");
                check1 = 0; check2 = 0;
                temp_cluster_no = icluster[0];
                assert (temp_cluster_no != 0);
                tempint = nclus[temp_cluster_no];
                                /* # of sequences grouped with QUERY */
                nclus[temp_cluster_no] = nclus[0];
                nclus[0] = tempint;
                for (s1 =0; s1 < nseq; s1++) {
                        if (icluster[s1] == temp_cluster_no) {
                                /* same cluster as seq 0*/
                                icluster[s1] = 0; check1++;
                        } else if (icluster[s1] == 0) {
                                /* in cluster 0, moving it to where
                                cluster number where query was */
                        icluster[s1] = temp_cluster_no;
                                check2++;
                        }
                }
                assert (check1 == nclus[0]); /* same # of sequences with query*/
                assert (check2 == nclus[temp_cluster_no]) ;
                        /* move all the sequences in cluster 0 to query */
        } /* end of if QUERY not infirst cluster */

	/*initialize cluster of seqs.  clusters_of_seqs is an array.
   each element in the array will contain the sequences in one cluster. */
      
   clusters_of_seqs =  calloc
				 (iclus, sizeof (struct ClusterOfSequences));

   for (i = 0; i < iclus; i++) {
	clusters_of_seqs[i].nseqs = 0;
   }
   for (s1 = 0; s1 < nseq; s1++)
   {
	cluster_num = icluster[s1];
	clusters_of_seqs[cluster_num].seqs[clusters_of_seqs[cluster_num].nseqs]
							 = seqs[s1];
        clusters_of_seqs[cluster_num].nseqs++;
   }
   for (i = 0 ; i < iclus ; i++) {
	assert (nclus[i] == clusters_of_seqs[i].nseqs);
   }
  *numclusters = iclus;
 printf ("Clumped %d sequences into %d clusters\n", nseq, iclus);
 /* keep this print because get # of clustersin awk statement */
  free(pairs);
printf ("exiting clump\n");
 return clusters_of_seqs;
}  /* end of cluster */

Block**
make_blocks_from_clumps (Sequence* seqs[MAXSEQ], /* IN */ 
			int nseqs, 		/* IN */
			double clumping_threshold, /* IN */
			 int* num_blocks)	/* OUT */
{
	struct ClusterOfSequences* cluster_of_seqs;
	Block ** array_of_blocks_made_from_clustered_seqs;
	int i;
	
	int num_clusters;

	cluster_of_seqs = cluster_seqs (clumping_threshold, seqs, nseqs,
				   &num_clusters);
	
	array_of_blocks_made_from_clustered_seqs =
		(Block**) calloc (num_clusters, sizeof (Block*));
	*num_blocks = num_clusters;

	 for (i = 0; i < num_clusters; i++) {
		array_of_blocks_made_from_clustered_seqs[i] =
		make_block (cluster_of_seqs[i].seqs[0]->length, /* width */
				 0, /* copy all of sequence */ 
				cluster_of_seqs[i].nseqs,
				cluster_of_seqs[i].seqs,
				 FALSE);
	}
	return array_of_blocks_made_from_clustered_seqs;		
} /* endof make_blocks_from_clumps */

/* returns an array of sequences taht contains consensus sequences 
doesn't work as of 2/26/01.   
*/
void
clump_into_consensus_seqs (Sequence* seqs[MAXSEQ], /* IN */
                        int nseqs,              /* IN */
                        double clumping_threshold, /* IN */
                         int* num_consensus_seqs, /* OUT*/
			Sequence* consensus_seqs[MAXSEQ])       /* OUT */
{
	int i;
        struct ClusterOfSequences* cluster_of_seqs;
        int num_clusters;
	Block* block; Matrix* pssm;

printf ("ok, an we even enter?\n");
        cluster_of_seqs = cluster_seqs (clumping_threshold, seqs, nseqs,
                                   &num_clusters);
	*num_consensus_seqs = num_clusters;
	printf ("number of clusters is %d\n", num_clusters);
	for (i = 0; i < num_clusters; i++) {
                /* save time, if there's only one sequence in the cluster then
  		obviously it's the consensus sequence */
		if (cluster_of_seqs[i].nseqs > 1) {
printf ("clusterlength %d\n", cluster_of_seqs[i].seqs[0]->length);
printf ("no of clusters %d\n", cluster_of_seqs[i].nseqs);
 
			block = make_block (cluster_of_seqs[i].seqs[0]->length, /* width */
                                 0, /* copy all of sequence */
                                cluster_of_seqs[i].nseqs,
                                cluster_of_seqs[i].seqs,
                                 FALSE);
        		pb_weights (block);
	        	pssm = PN_block_to_matrix (block, 2);
			consensus_seqs[i] = get_consensus (pssm);
			free_block(block); free_matrix (pssm);	
		} else {
			consensus_seqs[i] = (Sequence*) malloc (sizeof (Sequence));
			consensus_seqs[i]->sequence = (Residue *) calloc (cluster_of_seqs[i].seqs[0]->length, sizeof (Residue));
			copy_sequence (consensus_seqs[i], 
				cluster_of_seqs[i].seqs[0]);
		}
		sprintf (consensus_seqs[i]->name, "CONSENSUS%d", i);
	}

	free_ClusterOfSequences (cluster_of_seqs,num_clusters); 	 
}


int
print_consensus_seqs_aligned  (double clus, Sequence* seqs[MAXSEQ], 
				int nseq, FILE* outfp,
			        FILE* keyoutfp, int option)
{

   int iclus, npair, threshold, s1, s2, l1, l2, px, i, i1, i2, first, j;
   int check1, check2, temp_cluster_no, tempint;
   int nclus[MAXSEQ], icluster[MAXSEQ], minclus, oldclus, width;
   int match, no_of_clumps;
   struct cluster_pair *pairs;
  int alignment_width;
  int nseqs_cluster;
  Sequence* seqs_cluster[MAXSEQ]; Block* block; Matrix* pssm;
        Sequence* consensus_sequence;
        char comment[SMALL_BUFF_LENGTH];

   no_of_clumps = 0;
   npair = nseq*(nseq-1)/2;
   pairs = (struct cluster_pair *) malloc(npair * sizeof(struct cluster_pair));

   /*    Compute scores for all possible pairs of sequences            */
   for (s1=0; s1<nseq-1; s1++)                  /* col = 0, n-2     */
   {
      for (s2=s1+1; s2<nseq; s2++)      /* row = col+1, n-1 */
      {
         alignment_width = seqs[0]->length;
         match = 0;
         px = INDEX(nseq, s1, s2);
         pairs[px].score = 0;
         pairs[px].cluster = -1;
         for (i=0; i< seqs[0]->length; i++)
         {
                if ( seqs[s1]->sequence[i] == seqs[s2]->sequence[i]) {
                        if (seqs[s1]->sequence[i] == aa_atob['X'] ||
                              seqs[s1]->sequence[i]== aa_atob['-']) {
                                alignment_width--;
                        } else {
                                match++;
                        }
                } /* end of if letters match */
          } /* end of for width */
         pairs[px].score = ((double) match)/ ((double) alignment_width);
      }  /* end of s2 */
   }  /* end of s1 */

   /*-------Cluster if score exceeds threshold by scanning cols (s1) */
   for (s1=0; s1<nseq; s1++)
   {
      icluster[s1] = -1;                        /* clear out old values */
      nclus[s1] = 0;
   }
   iclus = 0;                                   /* cluster number */
   for (s1=0; s1<nseq-1; s1++)                  /* col = 0, n-2     */
      for (s2=s1+1; s2<nseq; s2++)      /* row = col+1, n-1 */
      {
         px = INDEX(nseq, s1, s2);
         if (pairs[px].score >= clus)      /*  cluster this pair */
         {
            if (icluster[s1] < 0)          /* s1 not yet clustered */
            {
               if (icluster[s2] < 0)       /* new cluster */
               {
                 icluster[s1] = iclus++;
                  icluster[s2] = icluster[s1];
                }
               else                             /* use s2's cluster  */
                  icluster[s1] =  icluster[s2];
            }
            /*  use s1's cluster if it has one and s2 doesn't */
            else if (icluster[s1] >= 0 && icluster[s2] < 0)
               icluster[s2] = icluster[s1];
            /* merge the two clusters into the lower number */
            else if (icluster[s1] >= 0 && icluster[s2] >= 0)
            {
               minclus = icluster[s1]; oldclus = icluster[s2];
               if (icluster[s2] < icluster[s1])
               {
                  minclus = icluster[s2]; oldclus = icluster[s1];
               }
               for (i1=0; i1<nseq; i1++)
                 if (icluster[i1] == oldclus)
                     icluster[i1] = minclus;
            }
         }  /* end of if pairs */
      }  /* end of s2 */

   /*---  Set ncluster, get rid of negative cluster numbers --*/
   for (s1=0; s1<nseq; s1++)
   {
      if (icluster[s1] < 0)
          icluster[s1] = iclus++;
   }
   for (s1=0; s1<nseq; s1++)
          nclus[icluster[s1]] += 1;

   i2 = 0;
/* printf ("no of clusters %d\n", iclus); */

   /* move QUERY to first cluster */
        if (icluster[0] != 0) { /* QUERY is not first in cluster */
printf ("entered if statment\n");
                check1 = 0; check2 = 0;
                temp_cluster_no = icluster[0];
                assert (temp_cluster_no != 0);
                tempint = nclus[temp_cluster_no];
                                /* # of sequences grouped with QUERY */
                nclus[temp_cluster_no] = nclus[0];
                nclus[0] = tempint;
                for (s1 =0; s1 < nseq; s1++) {
                        if (icluster[s1] == temp_cluster_no) {
                                /* same cluster as seq 0*/
                                icluster[s1] = 0; check1++;
                        } else if (icluster[s1] == 0) {
                                /* in cluster 0, moving it to where
                                cluster number where query was */
                        icluster[s1] = temp_cluster_no;
                                check2++;
                        }
                }
                assert (check1 == nclus[0]); /* same # of sequences with query*/
                assert (check2 == nclus[temp_cluster_no]) ;
                        /* move all the sequences in cluster 0 to query */
        } /* end of if QUERY not infirst cluster */
   for (i=0; i< iclus; i++)
   {
      if (nclus[i] > 1 || (nclus[i] > 0 && option == 0))
         /* print clumps that have 2 or more sequences or
        if option 0, print out clumps with 1 or more sequences */
      {
        i2++;
         nseqs_cluster = 0;
        for (s1=0; s1 < nseq; s1++)
         {
            if (icluster[s1] == i)
            {
                seqs_cluster[nseqs_cluster] = seqs[s1];
                nseqs_cluster++;
/*              output_sequence (seqs[s1], outfp); */
            }
         } /* end for nseqs */
        /* get consensus for cluster i */
        block = make_block (seqs[0]->length, 0, nseqs_cluster, seqs_cluster, FALSE);
        pb_weights (block);
        pssm = PN_block_to_matrix (block, 2);
        consensus_sequence = get_consensus (pssm);
/*              fprintf (outfp, "CONSENSUS for cluster %d \n", i); */
        /* August 28, 2000 added 1 to i below so that name starts as
        CONSENSUS1, CONSENSUS2, ... because QUERY will be represented
        by itself */
        sprintf (comment, "%d", i);
        strcat (consensus_sequence->name, comment);
        convert_gap_to_X (consensus_sequence);
        if (seq_not_all_Xes (consensus_sequence) ) {
                output_sequence (consensus_sequence, outfp);
                no_of_clumps++;
        }

/* now print out consensus key */
        fprintf (keyoutfp, "%s\n", consensus_sequence->name);
        for (j = 0; j < nseqs_cluster; j++) {
                fprintf (keyoutfp, "%s %s\n", seqs_cluster[j]->name, seqs_cluster[j]->info);
        }
        fprintf (keyoutfp, "\n");
        free_block (block);
        free_matrix (pssm);
        free_sequence (consensus_sequence);
     } /* end if nclus > 1 || nclus[i] > 0 && option == 0 */
   }


 printf ("Clumped %d sequences into %d clusters\n", nseq, i2);
 /* keep this print because get # of clustersin awk statement */
  free(pairs);
printf ("exiting clump\n");
 return no_of_clumps;
}  /* end of cluster */

#endif
