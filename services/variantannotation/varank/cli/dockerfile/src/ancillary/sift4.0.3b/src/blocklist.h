/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */   

#ifndef _BLOCKLIST_H
#define _BLOCKLIST_H

extern FILE* errorfp;

/*  first entry has no block, just nblock & totwidth, other entries
    in list have just pointers to the blocks */
struct block_listp {             /* list of blocks for a family */
   int nblock;                          /* number of blocks in family */
   int nseq;                            /* number of sequences in blocks */
   int totwidth;                        /* total width of the blocks */
   char block_num[10];
   Block *block;
   struct block_listp *next;
};
typedef struct block_listp Block_List;
 
Block_List *make_block_list(char *);
void insert_blist();
void free_blist();

#endif
