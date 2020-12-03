/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */ 

#ifndef _STRINGLIST_H
#define _STRINGLIST_H

#define KEY_WIDTH 50

extern FILE* errorfp;

struct string_listp {             /* list of blocks for a family */
   int no_of_elements;                  /* number of elements */
   char title[KEY_WIDTH];
   char key[KEY_WIDTH];
   struct string_listp *next;
};
typedef struct string_listp StringList;
 
StringList* make_string_list();
void insert_stringlist(StringList* list, char key[KEY_WIDTH]);
void free_stringlist();
void output_stringlist(StringList* list, FILE* outfp);

#endif
