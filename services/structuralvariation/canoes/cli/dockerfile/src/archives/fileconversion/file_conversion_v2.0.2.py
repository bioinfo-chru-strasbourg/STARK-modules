
##################################
###      Programme Python      ###
###     file_conversion.py     ###
###   Author : Céline Besnard  ###
##################################
"""This module does the conversion from one file format to another."""


####################################
### import of external functions ###
####################################
import argparse
import csv
import os
import re
import subprocess
import sys
import time
from os.path import join as osj
from Bio import SeqIO
import pandas as pd

##############################
## Release v2.0.1 07/10/2021: tsv2vcf Jean-Baptiste Lamouche Fix spaces not allowed in INFO field regex add, contig check in reference fix, sampleID as col name of GT etc 
## Release v2.0.2 29/10/2021: tsv2vcf version old AnnotSV Jean-Baptiste Lamouche Fix col name DUP/INS etc has to be SV type, generate edit tmp file in output folder, minors fixes

##############################
### Recovery of fasta file ###
##############################
def fasta():
    """Recovery of sequences from fasta file"""
    reference = open(REF, 'r')  # Open reference file
    print("Recovery of fasta file")
    seq = {}
    for record in SeqIO.parse(reference, "fasta"):
        seq[str(record.id)] = str(record.seq)
        # Dictionary : e.g key:chr1, value:sequence
    reference.close()
    return seq



##############################
### Calling a bash command ###
##############################
def get_stdout(cmd, debug=False):
    """Input is a string of a shell cmd returns a list of strings,
        where every string is a line from the stdout."""

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out = p.stdout.readlines()
    if debug:
        print(out)
    out_list = [l.rstrip() for l in out]
    err = p.stderr.readlines()
    if debug:
        print(err)
    err_list = [l.rstrip() for l in err]
    assert len(err_list) == 0, print("ERROR: len(err_list)>0. cmd=", cmd, "err_list=", err_list)
    return out_list


#################
### Sort file ###
#################
def sort_output(output):
    """Allows you to sort a file, input : file to sort"""
    print('Sort output file')
    tmp_name = output + "_unsorted"
    os.rename(output, tmp_name)
    get_stdout("grep '^#' " + tmp_name + " > " + output)
    cmd2 = "grep -v '^#' "+ tmp_name + " | sort -k1,1V -k2,2n >> " + output
    get_stdout(cmd2)
    #os.remove(tmp_name)



############################################
### Delete comment lines added by Canoes ###
############################################
def delete_line(file_temp) :
    """Delete comment lines added by Canoes, or delete lines in VCF's header"""
    edit = open(file_temp, 'w+')
    print("Delete header line")
    csv_writer = csv.writer(edit, delimiter = "\t")
    with open(FILE, 'r') as csvfile :
        csv_reader = csv.reader(csvfile, delimiter = '\t')
        for row in csv_reader:
            if re.search('##', str(row)):  # Line start with ## in VCF file is header
                continue
            csv_writer.writerow(row)
    edit.close()



###########################################
### Recovery of header line in bed file ###
###########################################
def recup_line(file_temp) :
    """Recovery lines in BED's header and create file without header"""
    edit = open(file_temp, 'w+')
    print("Delete header line")
    header = False
    csv_writer = csv.writer(edit, delimiter = "\t")
    with open(FILE, 'r') as csvfile :
        csv_reader = csv.reader(csvfile, delimiter = '\t')
        for row in csv_reader:
            if re.search(r'^#', str(row[0])) : # Line start with # in VCF file is header line
                csv_writer.writerow(row[1:])
                header = True
            elif re.search(r'^chr', str(row[0])) :
                csv_writer.writerow(row)
    edit.close()
    return header



###############################################
### Checking the minimum number of columns  ###
###############################################
def column_verification(dataf) :
    """Checking the minimum number of columns to be able to make the conversion"""
    if len(dataf.columns) < 5 :
        sys.exit("Convertion is not possible because informations missing")



########################################
### Retrieving useful variable names ###
########################################
def bed_name(dataf) :
    """Retrieving useful variable names, input is dataframe"""
    print("Retrieving useful variable names")
    for col in dataf.columns :
        if re.search(r".*chr.*", col, re.IGNORECASE) :
        # regular expression : Ignore the case
        # one or more characters (.+) followed by 'V', then any character (.),
        # then 'start' and this is end of string ($)
            chr = re.search(r".*chr.*", col, re.IGNORECASE).group(0)
        elif re.search(r".*start.*", col, re.IGNORECASE) :
            start = re.search(r".*start.*", col,re.IGNORECASE).group(0)
        elif re.search(r".*end.*",col, re.IGNORECASE) :
            end = re.search(r".*end.*",col, re.IGNORECASE).group(0)
        elif re.search(r".*genotype.*", col, re.IGNORECASE) :
            genotype = re.search(r".*genotype.*", col, re.IGNORECASE).group(0)
        elif re.search(r".*type.*", col, re.IGNORECASE) :
            type = re.search(r".*type.*", col, re.IGNORECASE).group(0)
    return chr, start, end, genotype, type


########################################
### Retrieving useful variable names ###
########################################
def name_of_var(dataf) :
    """Retrieving useful variable names, input is dataframe"""
    print("Retrieving useful variable names")
    for col in dataf.columns :
        if re.search(r".+V.start", col, re.IGNORECASE) :
        # regular expression : Ignore the case
        # one or more characters (.+) followed by 'V', then any character (.),
        # then 'start' and this is end of string ($)
            sv_start = re.search(r".+V.start$", col, re.IGNORECASE).group(0)
        elif re.search(r".+V.chrom", col, re.IGNORECASE) :
            sv_chrom = re.search(r".+V.chrom", col,re.IGNORECASE).group(0)
        elif re.search(r".+V.type",col, re.IGNORECASE) :
            sv_type = re.search(r"SV.type$",col, re.IGNORECASE).group(0)
        elif re.search(r".+V.end", col, re.IGNORECASE) :
            sv_end = re.search(r".+V.end$", col, re.IGNORECASE).group(0)
    return sv_type, sv_chrom, sv_start, sv_end



#########################
### sample info field ###
#########################
def sample_of_header(dataf, sv_type) :
    """Recovery of name of samples, input is dataframe and name of column sv_type"""
    print("Recovery of name of samples")
    print("#[INFO] df columns", dataf.columns)
    #for index, name in enumerate(dataf.columns):
    #    if name == sv_type :
    #        i = index + 1
    if 'Samples_ID' in dataf.columns and not dataf['Samples_ID'].isnull().values.any():
        i = list(dataf.columns).index('Samples_ID')
        sample_id = dataf[dataf.columns[i]].unique() # Composite list of samples (without duplicates)
    
        list_sample_id = sample_id.tolist()
        print("#[INFO] list_sample_id", list_sample_id)
        unique_sample = list(set((','.join(list_sample_id).split(','))),)  # Transformation into a list
        
    elif 'FORMAT' in dataf.columns:
        i = list(dataf.columns).index('FORMAT')+1
        unique_sample = (dataf.columns[i],)
        print("#[INFO] Sample Unique ",unique_sample)
    else:
        print("#[ERROR] Sample is empty")
        exit()
    print("#[INFO] index ", i)

    return i, unique_sample
    #print("#[INFO] col sampleID ", i)
    #print("#[INFO] SampleID normally", dataf['Samples_ID'])
    #print("#[INFO] Sample firs row", dataf.iloc[1, :])
    #print(dataf.head())
    



#########################
### Remove duplicates ###
#########################
def deduplicate(dataf, sv_type):
    """remove duplicates, input is dataframe and name of column sv_type"""
    index_sample, list_sample = sample_of_header(dataf, sv_type)
    print("Remove duplicates")
    print("Recovery of ID")
    dict = {}
    for column in dataf :
        if column == dataf.columns[index_sample] :
            datafdrop = dataf.drop(columns=[column])
            datafdrop.drop_duplicates(keep='first', inplace=True, ignore_index=True)
        if column == dataf.columns[0]:
            for sple in list_sample :
                print("#[INFO] SAMPLES CHECK", sple)
                loc = []
                # data_loc = dataf.loc[dataf[dataf.columns[index_sample]]==s, dataf.columns[0]]
                data_loc = dataf.loc[dataf[dataf.columns[index_sample]].str.contains(sple) == True, dataf.columns[0]]
                data_loc = data_loc.reset_index(drop=True)
                for i,j in enumerate(data_loc) :
                    loc.append(data_loc[i])
                dict[sple] = []
    print("#[INFO] DICT ",dict)
    return datafdrop, dict, list_sample



#########################
### header info field ###
#########################
def info_of_header(dataframe, sv_chrom, sv_start):
    """Obtaining information about the 'INFO' field of the header, input is dataframe"""
    info = ''
    dataf = dataframe.drop(columns=[sv_chrom, sv_start])
    if 'REF' in dataf.columns :
        dataf = dataframe.drop(columns=['REF', 'ALT'])
    for key,value in dico_info().items() :
        if key in dataf.columns :
            if str(dataf[key].dtypes) == 'object' :    #type of variables in the key column
                type = 'String'
            elif str(dataf[key].dtypes) == 'float64' :
                type = 'Float'
            else :
                type = 'Integer'
            info = info + '##INFO=<ID=' + key + ",Number=1,Type=" + type + ',Description="' + value + '">\n'
    info = info[:-1]      ##delete the last character of the string (\n)
    return info



########################
### header alt field ###
########################
def alt_of_header(dataframe, sv_type):
    """Obtaining information about the 'ALT' field of the header, input is dataframe"""
    alt = ''
    for key,value in dico_alt().items() :
        for d in dataframe[sv_type].unique() :
            if d.startswith('<') :
                d = d[1:-1]
            if key == d :
                alt = alt + '##ALT=<ID=' + key + ',Description="' + value + '">\n'
    alt = alt[:-1]      ##delete the last character of the string (\n)
    return alt



###########################
### header contig field ###
###########################
def contig_of_header(dataframe, sv_chrom, seq):
    """Obtaining information about the 'contig' field of the header, input is dataframe"""
    contig = ''
    chr = dataframe[sv_chrom].unique()
    chr = [ 'chr'+str(items) for items in chr]
    #print("[INFO] chr ", chr)
    for key, value in seq.items():
        #print("#[INFO] keys seq Reference ",key)
        if key in chr:
            contig = contig + '##contig=<ID=' + key + ',length=' + str(len(value)) + ',assembly=' + os.path.abspath(REF) +'>\n'
    #print("[INFO] Contig String ",contig )
    contig = contig[:-1]     # delete the last character of the string (\n)
    return contig



#############################################
### column description for the info field ###
#############################################
#def dico_info():
#    """DICTIONARY (column ID -> key, description -> value)"""
#    info_header = {}
#    with open("dico_head_" + version + ".txt",'r') as file :
#        line = file.readline().strip().decode('utf-8')
#        while line:
#            element = line.split(':')
#            info_header[element[0]] = element[1]
#            line = file.readline().strip()
#    return info_header

def dico_info():
    info_header = {}
    with open(osj(workingdirectory, "dico_head_" + version + ".txt"),'r') as file:
        for lines in file:
            #print("#[INFO] Lines without decode ", lines)
            lines = lines.strip().split(':')
            #print("#[INFO] Lines after decode ", lines)
            info_header[lines[0]] = lines[1]
    return info_header

############################################
### column description for the alt field ###
############################################
def dico_alt():
    """DICTIONARY (column ALT -> key, description -> value)"""
    #os.chdir("/STARK/input/tool/snakemake_scramble/scripts")
    alt_header = {}
    with open(osj(workingdirectory, "dico_alt.txt"),'r') as file :
        line = file.readline().strip()
        while line:
            element = line.split('=')
            alt_header[element[0]] = element[1]
            line = file.readline().strip()
    return alt_header



#########################
### creation of header ##
#########################
def header_generation(dataframe, sv_type, sv_chrom, seq) :
    """Creation of header in VCF file, input is dataframe and list of sample"""
    output = open(out,'w')
    print("Create header of VCF")
    date_of_the_day = time.strftime("%d/%m/%Y")        # dd/mm/yyyy format
    date = "##fileDate=%s" % date_of_the_day
    vcf_format = "##fileformat=VCFv4.3"       # VCF version
    source = "##source=VCF_file"
    input_file = "##InputFile=%s" % os.path.abspath(FILE)
    reference = "##reference=%s" % os.path.abspath(REF)        # reference genome
    filter = """##FILTER=<ID=PASS,Description="Passed filter">"""
    alt = alt_of_header(dataframe, sv_type)
    inf = """##INFO=<ID=SV_END,Number=1,Type=Integer,Description="Ending position of the SV in the chromosome">
##INFO=<ID=SV_TYPE,Number=1,Type=String,Description="Type of SV (DEL, DUP, ...)">
##INFO=<ID=SV_LEN,Number=1,Type=Integer,Description="Length of SV (bp)(deletions have negative values)">"""    # INFO declarations for header
    contig = contig_of_header(dataframe, sv_chrom, seq)
    format = """##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"""

    vcf_header_list = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'Sample']
    vcf_header = '\t'.join(vcf_header_list)
    output.write('\n'.join([vcf_format, date, source, input_file, reference, contig, alt, inf, filter, format, vcf_header]) + '\n')
    output.close()




#########################
### creation of header ##
#########################
def header(dataframe, unique_sample, sv_type, sv_chrom, sv_start, seq) :
    """Creation of header in VCF file, input is dataframe and list of sample"""
    output = open(out,'w')
    print("Create header of VCF")
    date_of_the_day = time.strftime("%d/%m/%Y")        # dd/mm/yyyy format
    date = "##fileDate=%s" % date_of_the_day
    vcf_format = "##fileformat=VCFv4.3"       # VCF version
    source = "##source=AnnotSV_annotation"
    input_file = "##InputFile=%s" % os.path.abspath(FILE)
    reference = "##reference=%s" % os.path.abspath(REF)        # reference genome
    filter = """##FILTER=<ID=PASS,Description="Passed filter">"""
    alt = alt_of_header(dataframe, sv_type)
    info = info_of_header(dataframe, sv_chrom, sv_start)    # INFO declarations for header
    contig = contig_of_header(dataframe, sv_chrom, seq)
    format = """##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"""

    vcf_header_list = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    vcf_header_return = vcf_header_list + list(unique_sample)
    vcf = '\t'.join(vcf_header_return)
    output.write('\n'.join([vcf_format, date, source, input_file, reference, contig, alt, info, format, vcf]) + '\n')

    output.close()
    return vcf_header_return



#######################
### creation of line ###
########################
def column_of_annotsv2vcf(dataf, colonne, dict, sv_type, sv_chrom, sv_start, seq) :
    """creation of output file, input are dataframe, output file, dictionnary of sample/IDs,
    name of column sv_type, sv_chrom, sv_start"""
    out = pd.DataFrame(columns = (colonne))
    print(out)
    print("Creation of column of VCF")
    for i in range(len(dataf)) :
        presence = {}
        ############ If the tsv was generated by a vcf ############
        if 'REF' in dataf.columns :
            chr = 'chr' + str(dataf[sv_chrom][i])
            ref = str(dataf['REF'][i])
            alt = str(dataf['ALT'][i])
            qual = str(dataf['QUAL'][i])
            id = str(dataf['ID'][i])
            pos = str(dataf[sv_start][i])
            datafdrop = dataf.drop(columns = [sv_chrom,sv_start, 'REF', 'ALT', 'FORMAT', 'FILTER', 'QUAL', 'ID'])    # Delete the columns already used
            for key, values in dict.items() :
                if key in datafdrop :
                    datafdrop = datafdrop.drop(columns = [key])
            datafdrop = datafdrop.drop(columns = ["INFO"])
        ############ else (the tsv was generated by a bed) ############
        else :
            chr = 'chr' + str(dataf[sv_chrom][i]) # e.g. chr1
            sequence = seq[chr]
            ref = sequence[int(dataf[sv_start][i])].upper()
            alt = '<' + str(dataf[sv_type][i]) + '>'
            qual = '.'
            id = '.'
            pos = str(dataf[sv_start][i] - 1) #O-based (on passe de VCF à BED)
            datafdrop = dataf.drop(columns = [sv_chrom, sv_start])     # Delete the columns already used
        datafdrop, informations, presence = column_info(datafdrop, dataf, i, presence, dict)
        informations = informations[1:]
        new_row = {'#CHROM':chr, 'POS':pos, 'ID':id, 'REF':ref, 'ALT':alt, 'QUAL':qual, 'FILTER':'PASS', 'INFO':informations, 'FORMAT':'GT'}
        for key, values in presence.items():
            if key in out.columns :
                new_row[key] = values
        out = out.append(new_row, ignore_index=True)
    print(out)
    out = out[1:]
    print("AFTER")
    print(out)

    return out



#####################################
### creation of info and genotype ###
#####################################
def column_info(datafdrop, dataf, i, presence, dict) :
    """creation of info field of output file, input are dataframe, dataframe without some columns,
     the iterator i, the dictionnary of presence and dictionnary of sample/IDs, NB: space are not allowed in INFO field"""
    informations = ''
    print("#[INFO] DATAF ", dataf)
    print("#[INFO] DATAFDROP", datafdrop.head())
    print("I value ", i)
    print("#[INFO]Presence ", presence)
    print("#[INFO] Dict ", dict)
    for column in datafdrop :
        print(column)
        if column == datafdrop.columns[0] :
            print("#[INFO] COLUMNS ", column)
            for key, values in dict.items() :
                #print("#[INFO] key, values ",key, values)
                #print(dataf.columns)
                if key in dataf.columns :
                    genotype = str(dataf[key][i]).split(':')
                    presence[key] = genotype[0]
                elif len(values) > 0:
                    if datafdrop[column][i] in values :
                        presence[key] = '1/1'
                else :
                    presence[key] = './.'
            informations = informations + ";" + column + '=' + str(datafdrop[column][i])
        elif str(datafdrop[column][i]) == 'nan' : # If empty, it is not considered
            continue
        elif str(datafdrop[column][i]).find(';') > -1 :
            datafdrop = datafdrop.replace(to_replace = datafdrop[column][i], value = ",".join(str(datafdrop[column][i]).split(';')))
            informations = informations + ";" + column + '=' + str(datafdrop[column][i])
        else :
            informations = informations + ";" + column + '=' + str(datafdrop[column][i])
 
    #remove space in informations field
    new_informations = re.sub(r',(\s{1})',",", informations)
    new_informations = re.sub(r'[0-9](\s{1})',"", new_informations)
    new_informations = re.sub(r'(\s{1})',"_", new_informations)
    return datafdrop, new_informations, presence



#######################
### creation of line ##
#######################
def column_of_annotsv(dataf, out, sv_type, sv_chrom, sv_start, sv_end):
    """creation of output file, input are dataframe, output file, dictionnary of sample/IDs,
    name of column sv_type, sv_chrom, sv_start, sv_end"""
    i = 0
    while i < len(dataf):
        for colonne in dataf.columns :
            # If the line is 'full' (annotation mode) then we write the row
            if dataf[colonne][i] == 'full' :
                chr = 'chr' + str(dataf[sv_chrom][i]) # e.g. chr1
                type = str(dataf[sv_type][i])
                if type.startswith('<') :
                    type = type[1:-1]
                if args.default == True :
                    new_row = {'chr' : chr, 'debut' : str(dataf[sv_start][i]), 'fin' : str(dataf[sv_end][i]), 'name' : chr + '_' + str(dataf[sv_start][i]) + '_' + str(dataf[sv_end][i]) + '_' + type, 'score' : '0', 'strand' : '+'}
                else :
                    new_row = {'chr' : chr, 'debut' : str(dataf[sv_start][i]), 'fin' : str(dataf[sv_end][i]), 'name' : chr + '_' + str(dataf[sv_start][i]) + '_' + str(dataf[sv_end][i]) + '_' + type}
                out = out.append(new_row, ignore_index = True) # Filling the dataframe
        i += 1
    return out



#######################
### creation of line ###
########################
def column_of_vcf(dataf):
    """creation of output file, input is dataframe"""
    if args.default == True :
        out = pd.DataFrame(columns = ['chr', 'debut', 'fin', 'name', 'score', 'strand']) # Création d'une dataframe
    else :
        out = pd.DataFrame(columns = ['chr', 'debut', 'fin', 'name']) # Création d'une dataframe
    for i in range(len(dataf)):
        chr = dataf["#CHROM"][i]
        start = dataf['POS'][i]
        end, type = split(dataf, i)
        if args.default == True :
            new_row = {'chr': chr, 'debut':start, 'fin':end, 'name':chr+'_'+str(start)+'_'+end+'_'+type, 'score':'0', 'strand' : '+'}
        else :
            new_row = {'chr': chr, 'debut':start, 'fin':end, 'name':chr+'_'+str(start)+'_'+end+'_'+type}
        out = out.append(new_row, ignore_index=True)
        i += 1
    return out



#######################
### creation of line ###
########################
def column_of_bed(dataf, seq, chr, start, end, genotype, type) :
    """creation of output file, input are dataframe, dictionnary of sequence,
    sv_chrom, sv_start, sv_end, gentotype and sv_type"""
    out = pd.DataFrame(columns=["#CHROM", 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'sample']) # Création d'une dataframe
    for i in range(len(dataf)):
        chromosome = dataf[chr][i]
        pos_start = dataf[start][i] + 1
        pos_end = dataf[end][i]
        type_var = dataf[type][i]
        geno = dataf[genotype][i]
        sequence = seq[chromosome]
        ref = sequence[int(pos_start)]
        length = int(pos_end) - int(pos_start)
        info = 'sv_type=' + type + ';sv_end=' + str(pos_end)+';SV_LEN='+str(length)
        new_row = {"#CHROM" : chromosome, 'POS' : pos_start, 'ID' : '.', 'REF' : ref, 'ALT' : type_var, 'QUAL':'.', 'FILTER' : '.', 'INFO' : info, 'FORMAT' : 'GT', 'sample': geno}
        out = out.append(new_row, ignore_index=True)
        i += 1
    print(out)
    return out



############################
### Recovery of variable ###
############################
def split(dataf, i):
    """Recovery of colonne in info field"""
    dict = {}
    info = dataf['INFO'][i]
    list_info = info.split(';')
    for l in list_info :
        if re.search(r".*end", l, re.IGNORECASE) :
            sv_end = re.search(r".*end", l, re.IGNORECASE).group(0)
        elif re.search(r".*type", l, re.IGNORECASE) :
            sv_type = re.search(r".*type", l, re.IGNORECASE).group(0)
        li = l.split('=')
        dict[li[0]] = li[1]
    if(dict[sv_end] == '' or dict[sv_type] == ''):
        sys.exit("sv_end or sv_type no found in 'INFO' column")
    return dict[sv_end], dict[sv_type]



###################################
### Verification of columns name ###
####################################
def column_name(dataf):
    """Verficiation of name of columns"""
    if '#CHROM' and 'POS' and 'INFO' in dataf.columns :
        pass
    else :
        sys.exit("It's not a real VCF")



##############
### tsv2bed ##
##############
def annotsv2bed() :
    """Conversion from tsv file (AnnotSV profile) to bed file"""
    print("#### AnnotSV2bed ####")
    delete_line('edit.tsv') # Call delete_line function
    print("Dataframe from AnnotSV data (tsv)")
    data = pd.read_csv('edit.tsv', sep = '\t')    # Dataframe creation
    dataf = pd.DataFrame(data)
    sv_type, sv_chrom, sv_start, sv_end = name_of_var(dataf) # Call name_of_var function
    dataf_2, dict, sample = deduplicate(dataf, sv_type) # Call deduplicate function
    if args.default == True :
        new_dataf = pd.DataFrame(columns = ['chr', 'debut', 'fin', 'name', 'score', 'strand']) # Création d'une dataframe
    else :
        new_dataf = pd.DataFrame(columns = ['chr', 'debut', 'fin', 'name']) # Création d'une dataframe
    dataf = column_of_annotsv(dataf_2, new_dataf, sv_type, sv_chrom, sv_start, sv_end)
    # Writing the dataframe to the output file
    dataf.to_csv(out, sep = '\t', index = False, header = False)
    os.remove('edit.tsv') # Deleting the temporary file
    sort_output(out)



##############
### vcf2bed ##
##############
def vcf2bed():
    """Conversion from vcf file to bed file"""
    print("#### vcf2bed ####")
    delete_line('edit.vcf')
    print("Dataframe from VCF data")
    # dataframe creation with file created in delete_line
    data = pd.read_csv('edit.vcf', sep = '\t')
    dataf = pd.DataFrame(data)
    column_name(dataf)
    dataf = column_of_vcf(dataf)
    # Writing the dataframe to the output file
    dataf.to_csv(out, sep = '\t', index = False, header = False)
    os.remove('edit.vcf')
    sort_output(out)



##############
### tsv2vcf ##
##############
def annotsv2vcf() :
    """Conversion from tsv file (AnnotSV profile) to vcf file"""
    print("#### AnnotSV2vcf ####")
    seq = fasta()
    delete_line(os.path.dirname(out)+'/edit.tsv')
    print("Dataframe from AnnotSV data (tsv)")
    data = pd.read_csv(os.path.dirname(out)+'/edit.tsv', sep = '\t')    #dataframe creation
    dataf = pd.DataFrame(data)
    print("#[INFO] dataf head", dataf.head())
    sv_type, sv_chrom, sv_start, sv_end = name_of_var(dataf)
    print("#[INFO] SV features: svtype: %s \n svchrom: %s \n svstart: %s \n sv_end: %s" % (sv_type, sv_chrom, sv_start, sv_end))
    dataf_2, dict, sample = deduplicate(dataf, sv_type)
    colonne = header(dataf_2,sample, sv_type, sv_chrom, sv_start, seq)
    dataf = column_of_annotsv2vcf(dataf_2, colonne, dict, sv_type, sv_chrom, sv_start, seq)
    dataf.to_csv(out, sep = '\t', index = False, header = False, mode='a', encoding="utf-8")
    os.remove(os.path.dirname(out)+'/edit.tsv')
    sort_output(out)



##############
### bed2vcf ##
##############
def bed2vcf () :
    """Conversion from tsv file (bed profile) to vcf file"""
    header = recup_line('edit.bed')
    print("Dataframe from VCF data")
    if header :  #dataframe creation with file created in delete_line
        data = pd.read_csv('edit.bed', sep = '\t', header = 0)
        dataf = pd.DataFrame(data)
    else :
        colonne = ['chr', 'start', 'end', 'type', 'genotype']
        if args.column :
            colnames = args.column.split(',')
            data = pd.read_csv('edit.bed',sep='\t', header=None)
            dataf = pd.DataFrame(data)
            j = 0
            for i in colnames :
                dataf.rename(columns = {int(i) : colonne[int(j)]}, inplace = True)
                j = j + 1
        else :
            sys.exit('No column names for bed2vcf, please look the README')
    column_verification(dataf)
    seq = fasta()
    sv_chr, sv_start, sv_end, genotype, sv_type = bed_name(dataf)
    header_generation(dataf, sv_type, sv_chr, seq)
    dataf = column_of_bed(dataf, seq, sv_chr, sv_start, sv_end, genotype, sv_type)
    # Writing the dataframe to the output file
    os.remove('edit.bed')
    dataf.to_csv(out, sep = '\t', mode = 'a', index = False, header = False)
    sort_output(out)


#####################
### fonction main ###
#####################
if __name__ == '__main__' :
    #################
    ### arguments ###
    #################
    parser = argparse.ArgumentParser(description = 'Transform file in other format')
    parser.add_argument('-c', '--convert', type = str, help = 'Type of conversion', choices = ['tsv2vcf', 'tsv2bed', 'vcf2bed', 'bed2vcf'] ,required = True)
    parser.add_argument('-p', '--profil', type = str, help = 'Profil of tsv', choices = ['AnnotSV', 'bed'])
    parser.add_argument('-i', '--infile', type = str, help = 'Tsv table of AnnotSV', required = True)
    parser.add_argument('-col', '--column', type = str, help = 'Number of column for bed2vcf : chr,start,stop,type,genotype')
    parser.add_argument('-d', '--default', type = bool, help = 'Strand and score by default', default = False)
    parser.add_argument('-ref', '--reference', type = str, help = 'Reference fasta file which table was build')
    parser.add_argument('-o', '--outfile', type = str, help = 'Vcf output file', required = True)
    parser.add_argument('-v', '--version', type = str, help = 'AnnotSV version', default = '3.0')
    parser.add_argument("-w","--workingdirectory", type = str, help = 'workingdirectory', default = '/app/src/fileconversion')
    args = parser.parse_args()
    workingdirectory = args.workingdirectory
    if os.path.exists(args.infile) :
        FILE = args.infile
    else :
        sys.exit("Input File doesn't exist")

    out = args.outfile

    convertion = args.convert

    defaut = args.default

    if convertion == 'tsv2vcf' :
        profil = args.profil
        if profil == 'AnnotSV' :
            if os.path.exists(args.reference) :
                REF = args.reference
            else :
                sys.exit("Fasta File doesn't exist")
            version = args.version
            annotsv2vcf()
        elif profil == 'bed' :
            if os.path.exists(args.reference) :
                REF = args.reference
                bed2vcf()
            else :
                sys.exit("Fasta File doesn't exist")
        else :
            sys.exit("Conversion isn't possible")

    elif convertion == 'tsv2bed' :
        profil = args.profil
        if profil == 'AnnotSV' :
            version = args.version
            annotsv2bed()
        elif profil == 'bed' :
            sys.exit('Convert bed to bed is not possible')

    elif convertion == 'vcf2bed' :
        vcf2bed()

    elif convertion == 'bed2vcf' :
        if os.path.exists(args.reference) :
            REF = args.reference
            bed2vcf()
        else :
            sys.exit("Fasta File doesn't exist")
