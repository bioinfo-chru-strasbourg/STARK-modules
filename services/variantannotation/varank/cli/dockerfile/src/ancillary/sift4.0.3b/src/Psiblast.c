/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */   

#ifndef _PSIBLAST_C_
#define _PSIBLAST_C_

#include "Psiblast.h"
#include "Alignment.c" /* for the function read_psiblast_header_until_first*/
#include "PN_convert.c"
#include "stringhash.c"

void process_sequence_line(Sequence* seq, char* buff);

Aligned_Pair*
read_psiblast_entry ( FILE* fp, char Buffer[LARGE_BUFF_LENGTH])
{
/* passing first line in entry that begins with >seq_id in Buffer */
        Aligned_Pair* alignment, *new_alignment;
	Aligned_Pair * cursor;
	
        alignment = initialize_Aligned_Pair (Buffer);
	cursor = alignment;

        while (fgets (Buffer,LARGE_BUFF_LENGTH, fp) &&
               ( strstr (Buffer, "Score") == NULL)) {
        /* keep reading until get score line */
        }

	reading_alignment_at_score_line (alignment, Buffer, fp);

	while (strstr (Buffer, "Score =") != NULL) {
		new_alignment =initialize_Aligned_Pair (Buffer);
		cursor->next = new_alignment;
		cursor = new_alignment;
	 	reading_alignment_at_score_line (cursor, Buffer, fp);		
	}
	return (alignment);

} /* end of read_psiblast entry*/

Aligned_Pair*
initialize_Aligned_Pair (char Buffer[LARGE_BUFF_LENGTH]) 
{
	Aligned_Pair* alignment;
	char* stringptr;
	int length;

       alignment = (Aligned_Pair*) malloc (sizeof (Aligned_Pair));

        assert (Buffer[0] == '>');
        stringptr = get_token (Buffer);
        length = strlen (stringptr);
        strncpy (alignment->subject->name, &Buffer[1], length-1);
        alignment->subject->name[length] = '\0';
        alignment->subject->type = AA_SEQ;
        alignment->query->type = AA_SEQ;
	alignment->next = NULL;

        strcpy (alignment->subject->info, &Buffer[length +1]);
	alignment->subject->length = 0;
	alignment->subject->max_length = 0;
	alignment->query->length = 0;
	alignment->query->max_length = 0;

	return alignment;

} /* end of initialize_Aligned_Pair*/

void
reading_alignment_at_score_line (Aligned_Pair* alignment, 
				char Buffer[LARGE_BUFF_LENGTH], FILE* fp)
{
        char* buff;
	int get_start_pos;

	buff = strstr(Buffer, "Score=");
        if (buff != NULL) {
                sscanf (buff, "Score = %d", alignment->score);
        } else {
                fprintf (errorfp, "Unable to read Score, parsing incorrect");
                exit (-1);
        }
        buff = strstr (Buffer, "Expect = ");
        if (buff != NULL) {
                sscanf (buff, "Expect = %lf", alignment->evalue);
        } else {
                fprintf (errorfp, "Unable to read e-value, parsing incorrect");
                exit (-1);
        }

        /* get Identities, Positives, Gaps line */
        fgets (Buffer, LARGE_BUFF_LENGTH, fp);
        fgets (Buffer, LARGE_BUFF_LENGTH, fp); /* get blank line */
        fgets (Buffer, LARGE_BUFF_LENGTH, fp); /* get first line */;
        get_start_pos = TRUE;
        while (Buffer != "\n") {
                read_4_alignment_lines (Buffer, alignment, fp, get_start_pos);
                get_start_pos = FALSE; /* read 1rst 4 alignment lines*/
                                        /* already got the starting pos*/
        }
        fgets (Buffer, LARGE_BUFF_LENGTH, fp);
        /* should be start of next alignment.  either beginning with > */
        /* or Score = if same protein, different alignment */
}

/*read 4_alignment_lines*/
/* pass in first line as Buffer, reads Query, subject alignment, and the
   following newline.  reads the next line into Buffer */
void
read_4_alignment_lines (char Buffer[LARGE_BUFF_LENGTH], Aligned_Pair* alignment,
                        FILE* fp, int get_start_pos)
{
        char* strptr;

        assert (strstr (Buffer, "Query:") != NULL);
        strptr = get_token (Buffer);
        assert (strcmp (strptr, "Query:") == 0);
        strptr = get_token (NULL);
        if (get_start_pos == TRUE) {
                alignment->query->position = atoi (strptr);
        }
        strptr = get_token (NULL); /* this is the alignment */
        process_sequence_line (alignment->query, strptr);
        fgets (Buffer, LARGE_BUFF_LENGTH, fp);
        fgets (Buffer, LARGE_BUFF_LENGTH, fp);
        strptr = get_token (Buffer);
        assert (strcmp (strptr, "Sbject:") == 0);
        strptr = get_token (NULL);
        if (get_start_pos == TRUE) {
                alignment->subject->position = atoi (strptr);
        }
        strptr = get_token (NULL); /* subject alignment */
        process_sequence_line (alignment->subject, strptr);
        fgets (Buffer, LARGE_BUFF_LENGTH, fp); /* read in newline */
        fgets (Buffer, LARGE_BUFF_LENGTH, fp); /* read in next line */

} /* end of read_4_alignment_lines */

/*
 * process_sequence_db_line
 *   Passed one line of the sequence portion of the database.  Reads the
 *   characters and converts them to the apropriate internal representation.
 *   Removes the excess spaces in the middle and ends.
 *   NOTE: this function adds the residues onto the end of the sequence
 *         passed to it.
 *   Parameters:
 *     Sequence *seq: the sequence to put the residues into.
 *     char *buff:    the string with the characters to residues to
 *                      add to the sequence
 *   Error codes: none
 */

void process_sequence_line(Sequence* seq, char* buff)
{
  int start_pos;
  int estimated_length;
  int saved_length;
  int num_entered;

  estimated_length = strlen(buff);
  saved_length = seq->length;

  /* if there are more residues than there is space for, resize. */
  /* Note: The buff may have characters that are not residues, this is */
  /*       ok, there will just be a little extra space */
  while ((saved_length + estimated_length) > seq->max_length) {
    resize_sequence(seq); 
  }

  /* temporarily increase the room to put the sequence into */
  start_pos = seq->length;
  seq->length = seq->max_length;

  num_entered = read_sequence(seq, seq->type, start_pos, buff);

  seq->length = saved_length + num_entered;

}

void
free_aligned_pair (Aligned_Pair* alignment)
{
	Aligned_Pair * cursor;

	cursor = alignment;

	while (cursor->next != NULL) {

		free_sequence (alignment->query);
		free_sequence (alignment->subject);
		free (alignment);
	}
	alignment = NULL;
	
}

void block_to_checkpoint_file (Block* block, FILE* fp, FILE* outfp) 
{

	Matrix* pssm;
	long ltemp;
	int pos;
	char ctemp[12];
	int aa; double dtemp;
printf ("in block to checkpoint file\n");
	pssm = PN_block_to_matrix (block, 20);
/*	print_matrix (pssm);*/ 
printf ("after conversion\n");
	ltemp = (long) block->width;
	fwrite (&ltemp, sizeof (long), 1, fp);
	for (pos = 0; pos < block->width; pos++)
	{
		ctemp[0] = toupper (aa_btoa[block->sequences[0].sequence[pos]]);
		fwrite(&ctemp, sizeof (char), 1, fp);
	}
	for (pos = 0; pos < block->width; pos++)
	{
		for (aa = 1; aa < 21; aa++) {
			dtemp = pssm->weights[aa][pos];
			fwrite (&dtemp, sizeof (double), 1, fp);
		}
	}
	free_matrix (pssm);
printf ("finished converting to checkpointsuccessful\n");
	
} /* end block_to_checkpoint_file */

void
psiblast_system_call (char chkpoint_filename[], char database[],
			 char result_file[], char query_seq_file[], FILE* outfp)
{
	char command_line[LARGE_BUFF_LENGTH * 4];
        char* ncbi_dir;

	ncbi_dir = getenv ("NCBI");

	sprintf (command_line, "%s/blastpgp -d %s -i %s -o %s -m 6 -R %s\n",
		ncbi_dir, database, query_seq_file, result_file, 
		chkpoint_filename);
	system (command_line);

}

void
rm_file (char filename[])
{
	char command_line[LARGE_BUFF_LENGTH];

	sprintf (command_line, "rm -f %s\n", filename);
	system ("unalias rm");
	system (command_line);
}

void 
formatdb_system_call (char database[])
{
	char command_line[LARGE_BUFF_LENGTH *2];
	char* ncbi_dir;

	ncbi_dir = getenv ("NCBI");
	sprintf (command_line, "%s/formatdb -i %s -o T -p T\n", ncbi_dir, database);
	printf ("formatting database command_line %s\n", command_line);
	system (command_line);
	printf ("finished formatting database\n");
} /* end formatdb_system_call */

/* opens psiblast result file and reads the top-scoring matches.  The top match
that isn't already in the sequence name hash is returned in next_seq_to_add.
*/
void
get_top_seq (char next_seq_to_add[], HashTable namehash, 
		char psiblastres_file[])
{
	FILE* psiblastfp;
	char name[KEY_WIDTH];
	char *strptr;
	char line[LARGE_BUFF_LENGTH];


	if ((psiblastfp = fopen (psiblastres_file, "r")) == NULL)
	{
		fprintf (errorfp, "can't open psiblast results %s\n", psiblastres_file);
		exit (-1);
	}


	if (read_psiblast_header_until_first_no_error (psiblastfp, FALSE) == -1) {
		next_seq_to_add[0] = '\0'; 
	} else {

	fgets (line, LINE_LEN, psiblastfp);
	strptr = strtok(line, " \t");
	/*strptr[strlen (strptr)] = '\0'; */
	assert ( strlen(strptr) < KEY_WIDTH);
	strncpy (name, strptr, KEY_WIDTH);
	name[KEY_WIDTH-1] = '\0';
	} 
	while (Exists(name, namehash)) {
		fgets (line, LINE_LEN, psiblastfp); 
		strptr = strtok(line, " \t\n");
        	if (strptr == NULL) {
			printf ("entered no line\n");
			name[0] = '\0';
		} else {
			strncpy (name, strptr, KEY_WIDTH);
			/* 	strptr[strlen(strptr)] = '\0' ; */
			name[KEY_WIDTH-1] = '\0';      
		}
	}
	strncpy (next_seq_to_add, name, KEY_WIDTH);
	next_seq_to_add[KEY_WIDTH -1] = '\0';

	fclose (psiblastfp);
 
} /* end get_top_seq */

#endif

