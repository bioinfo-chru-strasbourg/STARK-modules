/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */   

#ifndef PSIBLAST_H_
#define PSIBLAST_H_

#include "blocksprogs.h"

#include <string.h>
extern FILE* errorfp;

Matrix* read_psiblast_pssm (FILE* psiblastfp, int length);

void read_psiblast_matrix_body (FILE* fp, Matrix* matrix);

int read_psiblast_header (FILE* fp) ;

char get_index_PSIBLAST_order (int index);


#endif

