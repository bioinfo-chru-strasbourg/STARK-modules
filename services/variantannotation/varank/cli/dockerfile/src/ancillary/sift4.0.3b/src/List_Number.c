/* COPYRIGHT 2000 , Fred Hutchinson Cancer Research Center */ 

#ifndef _LIST_NUMBER_C_
#define _LIST_NUMBER_C_

#include "List_Number.h"

NumNodePtr make_num_list()
{
	NumNodePtr new;

	new = (NumNodePtr) malloc (sizeof (struct NumNode));
	new->beg = 0;
	new->end = 0;
	new->pos = 0;

	new->no_of_elements = 0;
	new->next = NULL;

	return new;
}

void insert_beg_end (NumNodePtr num_list, int beg, int end)
{
	NumNodePtr cursor;
	NumNodePtr new;

	cursor = num_list;
	while (cursor->next != NULL) {
		cursor = cursor->next;
	}
	new = make_num_list();
	new->beg = beg;
	new->end = end;
	cursor->next = new;
	num_list->no_of_elements++;
} /* end insert_beg_end */

void free_NumNode (NumNodePtr node)
{
	NumNodePtr cursor;
	NumNodePtr temp;

	assert (node != NULL);

	cursor = node->next;
	node->next = NULL;
	while (cursor != NULL) {
		temp = cursor->next;
		free (cursor);
		cursor = temp;
	}
	free (node);
	node = NULL;

} /* end free_NumNode */

NumNodePtr
read_region_file (FILE* fp)
{
	char line[LARGE_BUFF_LENGTH];
	int beg, end;
	NumNodePtr list;
	char* ptr;

	list = make_num_list();
	while (!feof (fp) && (fgets (line, LARGE_BUFF_LENGTH, fp) != NULL) ) {
		ptr = strtok (line, " \t\n\r\0");
		beg = atoi (ptr);
		ptr = strtok (NULL, " \t\n\r\0");
		end = atoi (ptr);
		insert_beg_end (list, beg, end);
	}
 	return list; 
} /* end read_region_file */

#endif
