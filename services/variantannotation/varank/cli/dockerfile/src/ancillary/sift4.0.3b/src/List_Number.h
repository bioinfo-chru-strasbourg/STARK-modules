/* COPYRIGHT 2000 , Fred Hutchinson Cancer Research Center */ 

#ifndef _LIST_NUMBER_H_
#define _LIST_NUMBER_H_

struct NumNode {
	int beg; 	/* things I might want to store */
	int end;
	int pos;
	
	int no_of_elements;
	struct NumNode* next;
};

typedef struct NumNode* NumNodePtr;

NumNodePtr make_num_list();
void insert_beg_end (NumNodePtr num_list, int beg, int end);
void free_NumNode (NumNodePtr node);

#endif
