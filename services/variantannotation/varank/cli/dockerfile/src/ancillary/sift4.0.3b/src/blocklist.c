/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */   

#ifndef BLOCKLIST_C_
#define BLOCKLIST_C_

#include "stdlib.h"
#include <string.h>
#include "blocklist.h"

#define MAXSEQ 400

/*======================================================================
routines for a list of blocks
========================================================================*/
Block_List *make_block_list(char* block_name)
{
   Block_List *new;

   new = (Block_List *) malloc (sizeof(Block_List));
   new->nblock = new->nseq = new->totwidth = 0;
   new->block = NULL;
   new->next = NULL;
   strcpy (new->block_num, block_name);
   new->block_num[7] = '\0';
   return(new);
}  /* end of make_blist */

void insert_blist(blist, block, desc)
Block_List *blist;
Block *block;
char* desc;
{
   char ctemp[SMALL_BUFF_LENGTH];
   Block_List *cur;
  
   assert (block != NULL);
   cur = blist;
   while (cur->next != NULL)
      cur = cur->next;
   strcpy (block->number, desc);
   cur->next = make_block_list(block->number);
   cur->next->block = block;
   blist->nblock += 1;
   if (block->num_sequences > blist->nseq)
      blist->nseq = block->num_sequences;
   blist->totwidth += block->sequences[0].length;

   /*------ Fill in some of the block information -------*/
   strcpy (block->id, desc);
   strcpy(block->de, desc);
   strcpy(block->motif, "UNK");
   sprintf(block->bl, "UNK motif; width=%d; seqs=%d;",
          block->width, block->num_sequences);
   strcpy (ctemp, desc);
   if      (blist->nblock ==  1) strcat(ctemp, "A");
   else if (blist->nblock ==  2) strcat(ctemp, "B");
   else if (blist->nblock ==  3) strcat(ctemp, "C");
   else if (blist->nblock ==  4) strcat(ctemp, "D");
   else if (blist->nblock ==  5) strcat(ctemp, "E");
   else if (blist->nblock ==  6) strcat(ctemp, "F");
   else if (blist->nblock ==  7) strcat(ctemp, "G");
   else if (blist->nblock ==  8) strcat(ctemp, "H");
   else if (blist->nblock ==  9) strcat(ctemp, "I");
   else if (blist->nblock == 10) strcat(ctemp, "J");
   else if (blist->nblock == 11) strcat(ctemp, "K");
   else if (blist->nblock == 12) strcat(ctemp, "L");
   else if (blist->nblock == 13) strcat(ctemp, "M");
   else if (blist->nblock == 14) strcat(ctemp, "N");
   else if (blist->nblock == 15) strcat(ctemp, "O");
   else if (blist->nblock == 16) strcat(ctemp, "P");
   else if (blist->nblock == 17) strcat(ctemp, "Q");
   else if (blist->nblock == 18) strcat(ctemp, "R");
   else if (blist->nblock == 19) strcat(ctemp, "S");
   else if (blist->nblock == 20) strcat(ctemp, "T");
   else if (blist->nblock == 21) strcat(ctemp, "U");
   else if (blist->nblock == 22) strcat(ctemp, "V");
   else if (blist->nblock == 23) strcat(ctemp, "W");
   else if (blist->nblock == 24) strcat(ctemp, "X");
   else if (blist->nblock == 25) strcat(ctemp, "Y");
   else if (blist->nblock == 26) strcat(ctemp, "Z");
   else strcat(ctemp, "*");

   sprintf(block->ac, "%s; distance from previous blocks=(0,0)",
					ctemp);

}  /* end of insert_blist */

void free_blist(blist)
Block_List *blist;
{
   Block_List *cur, *last;
   Block* tempblock;
   Block_List* P, *Tmp;
 
   assert ( blist != NULL);
 
   P = blist->next; /* Header assumed */
   blist->next = NULL;
   while (P != NULL) {
        Tmp = P->next;
        free_block (P->block);
        free(P);
        P = Tmp;
   }
   free(blist);
   blist = NULL;
    
 
}  /* end of free_blist */


#endif

