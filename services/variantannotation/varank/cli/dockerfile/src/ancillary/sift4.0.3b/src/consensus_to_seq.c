/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */ 

/* consensus_to_seq.c */
/* takes in sequence list and key list (key list has sequence names and
corresponding sequences with it) and prints out the corresponding sequences

Ex.  sequence list has CONSENSUS 0, CONSENSUS 3, CONSENSUS 37, key list
has
CONSENSUS0
QUERY
HBB_HUMAN
Q14473
...
CONSENSUS3
Q03901
HBD_GALCR

and will print out

QUERY
HBB_HUMAN
Q14473
Q03901
HBD_GALCR

*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "stringlist.c"

#define MAX_SEQ 400
#define LINE_LEN 400
#define true 1
#define false 0

void getargs (int argc, char* argv[], FILE** seqfp, FILE** keyfp, FILE** outfp);


/* MAIN */

int
main (int argc, char* argv[]) 

{

	StringList* consensus[MAX_SEQ]; /* array of lists of strings , with
					the title being consensus0, etc.
				and elements in the list being the items*/
	int i, no_of_items, done;
	FILE* keyfp, *seqfp, *outfp;
	char line[LINE_LEN];
	char* strptr;
	char temp_key[KEY_WIDTH];

	no_of_items = 0;

	getargs (argc, argv, &seqfp, &keyfp, &outfp);	

	while (!feof (keyfp)  && fgets (line, LINE_LEN, keyfp) != NULL) {
		strptr = strtok (line, " \t\r\n");
		strptr[strlen(strptr)] = '\0';
		consensus[no_of_items] = make_string_list();
		strcpy (consensus[no_of_items]->title, strptr);
		while (!feof (keyfp) && 
			fgets (line, LINE_LEN, keyfp) && line[0] != '\n') {
			strptr = strtok (line, " \t\n\r");
			strptr[strlen(strptr)] = '\0'; 
			strcpy (temp_key, strptr);
			insert_stringlist(consensus[no_of_items], temp_key); 
/*			insert_stringlist (consensus[no_of_items], line); */
		}
		/* done with this consensus and its corresponding elements*/
		/* go on to next list */
		no_of_items++;
	}
	fclose (keyfp); 
	while (!feof (seqfp) && fgets(line, LINE_LEN, seqfp) != NULL) {
		i = 0;
                strptr = strtok (line, " \t\r\n");
		strptr[strlen(strptr)] = '\0'; done = false;
		while ( i < no_of_items && !done) {
		if (strcmp (consensus[i]->title, strptr) == 0) {
				done = true; i--;
			}
			i++;
		}
		if (i == no_of_items) {
			printf ("couldn't find %s\n", strptr);
			exit (-1);
		}
		output_stringlist (consensus[i], outfp);
	}
	fclose (seqfp); fclose (outfp);
	for (i = 0; i < no_of_items; i++) {
		free_stringlist (consensus[i]);
	}
	exit (0);
}


void getargs
(int argc, char* argv[], FILE** seqfp, FILE** keyfp, FILE** outfp)
{
        char seqfilename[LINE_LEN];
        char outfilename[LINE_LEN];
        char keyfilename[LINE_LEN];

	if (argc < 3) {
		printf ("consensus_to_seqs.c: Looks up elements corresonding ");
		printf ("to key and prints it out\n");
	}

	if (argc > 1) strcpy (seqfilename, argv[1]);
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

	if (argc > 2) strcpy (keyfilename, argv[2]) ;
	else
	{
		printf ("Enter file with keys and elements\n");
		gets (keyfilename);
	}
	if ((*keyfp = fopen (keyfilename, "r")) == NULL) 
	{
		printf ("cannot open %s\n", keyfilename);
		exit (-1);
	}
	if (argc > 3) strcpy (outfilename, argv[3]);
	else
	{
		printf ("Enter name of outfile\n");
		gets (outfilename);
	}
	if ( (*outfp = fopen (outfilename, "w")) == NULL)
	{
		printf ("Cannot open %s\n", outfilename);
		exit (-1);
	}

} /* end of getargs */
	 
