/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */ 

/* stringhash.h
*/


        #ifndef _STRINGHASH_H_
        #define _STRINGHASH_H_

#define KEY_WIDTH 30

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern FILE* errorfp;

        typedef unsigned int Index;

	struct ListNode;
	typedef struct ListNode *Position;
        
	struct HashTbl;
        typedef struct HashTbl *HashTable;

        HashTable InitializeTable( int TableSize );
      
	void DestroyTable( HashTable H );
        
	Position Find( char* Key, HashTable H );
	
	int Exists (char * Key, HashTable H);
        
	void Insert( char* Key, HashTable H );
       

        #endif  /* _Stringhash_H */

