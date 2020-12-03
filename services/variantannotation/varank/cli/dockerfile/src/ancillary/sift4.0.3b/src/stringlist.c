/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */ 

#ifndef STRINGLIST_C_
#define STRINGLIST_C_

#include "stdlib.h"
#include <string.h>
#include "stringlist.h"


/*======================================================================
routines for a list of strings
========================================================================*/
StringList* make_string_list()
{
   StringList* new;

   new = (StringList*) malloc (sizeof(StringList));
   new->no_of_elements = 0; 
   new->next = NULL;
   return(new);
}  /* end of make_blist */

void insert_stringlist(StringList* list, char key[KEY_WIDTH])
{
   StringList *cur;
  
   cur = list;
   while (cur->next != NULL)
      cur = cur->next;
   cur->next = make_string_list();
   strcpy (cur->next->key, key);
   list->no_of_elements++;

}  /* end of insert_stringlist */

void free_stringlist(StringList* list)
{
   StringList *cur, *last;
   StringList* P, *Tmp;
 
   assert ( list != NULL);
 
   P = list->next; /* Header assumed */
   list->next = NULL;
   while (P != NULL) {
        Tmp = P->next;
        free(P);
        P = Tmp;
   }
   free(list);
   list = NULL;
    
 
}  /* end of free_blist */

void output_stringlist (StringList* list, FILE* outfp) 
{
	StringList* cur;

	cur = list->next;
	assert (cur != NULL); /* there has to be at least one item in the list*/
	while (cur != NULL) {
		fprintf (outfp, "%s\n", cur->key);
		cur = cur->next;
	}
}

#endif

