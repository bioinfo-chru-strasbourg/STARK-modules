#!/bin/bash
## STARK application GENOME

# DEFAULT ENV
#############
# Application heritage (default if no heritage)
#source_app $CONFIG_DEFAULT_APP


# DEJAVU PARAMETERS
#####################
# Spécific variables for DEJAVU script
# Don't forget to use "export"

# DOCKER DATABASES SERVICE DEJAVU listener SLEEP (1day=24hours=86400)
export SLEEP=86400


# STARK PARAMETERS
####################
# No need to use export for variables of STARK application (processed within script)

# All parameters as STARK application configuration file can be changed
# Here are some useful


# THREADS (default AUTO)
# Number of threads to use for the analysis
# AUTO will considere CORE-1 threads to use
# The number of threads need to be between 1 and the total number of cores available (autoadjusting if bad value, default AUTO i.e. all threads)
#THREADS=AUTO


# ANNOTATION
# Default annotation with HOWARD for intermediate VCF (for each caller) used by default with annotation rule "howard"
#ANNOTATION_TYPE="core,frequency,score,annotation,prediction,snpeff,snpeff_hgvs" "core,symbol,location,outcome,hgvs,snpeff,snpeff_hgvs,snpeff_split"
# Default annotation with HOWARD for report
HOWARD_ANNOTATION_REPORT="core,frequency,score,annotation,prediction,snpeff,snpeff_hgvs,snpeff_split"


# CALCULATION
# Default calculation with HOWARD for final VCF report
#HOWARD_CALCULATION_REPORT="FindByPipelines,GenotypeConcordance,VAF,VAF_STATS,DP_STATS,VARTYPE,NOMEN,BARCODE"
# List of annotation fields to extract NOMEN annotation (default 'hgvs', see HOWARD docs)
HOWARD_NOMEN_FIELDS="snpeff_hgvs"
