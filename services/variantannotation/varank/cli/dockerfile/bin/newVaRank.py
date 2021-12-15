# !/usr/bin/python
# !/bin/sh
# !/usr/tcl/bin
# -*- coding: utf-8 -*-
###############################
##							 ##
##		VARANK ANALYSIS		 ##
## need alamut-batch license ##
##	Author : Mateusz RAUCH	 ##
##							 ##
###############################

import os
import argparse
import sys
import subprocess
import shutil
import glob
from datetime import datetime
import re
import pandas


def main():
	launcher()


def launcher():
	original_umask = os.umask(0o000)

	args = parse_args()
	args_checker(args)
	varank_config_folder_checker(args)
	environment_checker(args)

	if not args.configure:
		varank_config_folder_checker(args)
		if args.run:
			launch_run_analysis(args)
		elif args.varank_folder:
			launch_folder_analysis(args)
		elif args.allsync18:
			pattern = pattern_generator(args)
			launch_universal_sync18(pattern)

	os.umask(original_umask)


def launch_universal_sync18(pattern):
	varank_processing_vcf_files = glob.glob(os.path.join(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_SERVICES'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_MODULE_NAME'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_SUBMODULE_NAME'], 'Archives', '*', '*', 'VCF', '*', '*.final.vcf*'))

	stark_archives_all_vcf_files = glob.glob(os.path.join(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_ARCHIVES'], '*', '*', '*', '*/*.final.vcf.gz'))
	for stark_archives_all_vcf_file in stark_archives_all_vcf_files:
		platform = stark_archives_all_vcf_file.split('/')[4]
		application = stark_archives_all_vcf_file.split('/')[5]
		run = stark_archives_all_vcf_file.split('/')[6]
		varank_processing_run_vcf_folder = os.path.join(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_SERVICES'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_MODULE_NAME'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_SUBMODULE_NAME'], 'Archives', platform, application, 'VCF', run)
		if not os.path.isdir(varank_processing_run_vcf_folder):
			os.makedirs(varank_processing_run_vcf_folder)
		subprocess.run(['rsync', '-rp', stark_archives_all_vcf_file, varank_processing_run_vcf_folder])

	for all_vcf in reversed(pattern):
		if all_vcf == '*/*.final.vcf.gz':
			pattern.remove(all_vcf)
		else:
			all_other_vcf_files = glob.glob(os.path.join(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_ARCHIVES'], '*', '*', '*', all_vcf))
			for all_other_vcf_file in all_other_vcf_files:
				platform = all_other_vcf_file.split('/')[4]
				application = all_other_vcf_file.split('/')[5]
				run = all_other_vcf_file.split('/')[6]
				varank_processing_run_vcf_folder = os.path.join(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_SERVICES'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_MODULE_NAME'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_SUBMODULE_NAME'], platform, application, 'VCF', run)
				if not os.path.isdir(varank_processing_run_vcf_folder):
					os.makedirs(varank_processing_run_vcf_folder)
				subprocess.run(['rsync', '-rp', all_other_vcf_file, varank_processing_run_vcf_folder])
			pattern.remove(all_vcf)

		for varank_processing_vcf_file in varank_processing_vcf_files:
			os.chmod(varank_processing_vcf_file, 0o777)


def vcf_synchronizer(pattern, run_repository, varank_processing_run_vcf_folder, varank_processing_folder):
	stark_vcf_files = glob.glob(os.path.join(run_repository, pattern[0]))
	for stark_vcf_file in stark_vcf_files:
		if not os.path.isdir(varank_processing_run_vcf_folder):
			os.makedirs(varank_processing_run_vcf_folder)
		subprocess.run(['rsync', '-rp', stark_vcf_file, varank_processing_run_vcf_folder])

	for vcf in reversed(pattern):
		if vcf == pattern[0]:
			pattern.remove(vcf)
		else:
			other_vcf_files = glob.glob(os.path.join(run_repository, vcf))
			for other_vcf_file in other_vcf_files:
				if not os.path.isdir(varank_processing_run_vcf_folder):
					os.makedirs(varank_processing_run_vcf_folder)
				subprocess.run(['rsync', '-rp', other_vcf_file, varank_processing_run_vcf_folder])
			pattern.remove(vcf)

	varank_processing_vcf_files = glob.glob(os.path.join(varank_processing_folder, 'VCF', '*', '*.final.vcf*'))
	for varank_processing_vcf_file in varank_processing_vcf_files:
		os.chmod(varank_processing_vcf_file, 0o777)


def varank_processing_folder_initializer(varank_processing_tsv_folder, varank_processing_folder):
	varank_processing_vcf_files_run_vcf_folder = glob.glob(os.path.join(varank_processing_folder, 'VCF', '*', '*.final.vcf*'))
	varank_processing_run_vcf_folder = glob.glob(os.path.join(varank_processing_folder, 'VCF', '*'))
	varank_processing_vcf_files_vcf_folder = glob.glob(os.path.join(varank_processing_folder, 'VCF', '*.final.vcf*'))
	varank_processing_vcf_files_main_folder = glob.glob(os.path.join(varank_processing_folder, '*.final.vcf*'))
	other_varank_processing_vcf_files = varank_processing_vcf_files_main_folder + varank_processing_vcf_files_vcf_folder

	deleted_samples_file = os.path.join(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_CONFIG'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_MODULE_NAME'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_SUBMODULE_NAME'], 'configfiles', 'deleted_samples_config.tsv')

	if os.path.isfile(deleted_samples_file):
		with open(deleted_samples_file, 'r') as read_file:
			next(read_file)
			for line in read_file.readlines():
				line = line.split('\t')
				if line[4] == os.path.basename(varank_processing_folder):
					vcf_file_to_delete = os.path.join(varank_processing_folder, 'VCF', line[1], line[0] + '.final.vcf.gz')
					if os.path.isfile(vcf_file_to_delete):
						os.remove(vcf_file_to_delete)
						print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: Deleted sample ' + line[0] + ' from run ' + line[1] + ' for this VaRank analysis')

	if len(varank_processing_vcf_files_run_vcf_folder) != 0:
		varank_processing_run_vcf_folder = sorted(varank_processing_run_vcf_folder)
		for varank_processing_run in varank_processing_run_vcf_folder:
			varank_processing_vcf_file_list = glob.glob(os.path.join(varank_processing_run, '*.final.vcf*'))
			for varank_processing_vcf_file in varank_processing_vcf_file_list:
				subprocess.run(['rsync', '-rp', varank_processing_vcf_file, varank_processing_folder])

	elif len(other_varank_processing_vcf_files) != 0:
		if not os.path.isdir(os.path.join(varank_processing_folder, 'VCF')):
			os.mkdir(os.path.join(varank_processing_folder, 'VCF'), 0o777)
		if not os.path.isdir(os.path.join(varank_processing_folder, 'VCF', 'DEFAULT')):
			os.mkdir(os.path.join(varank_processing_folder, 'VCF', 'DEFAULT'), 0o777)
		varank_processing_vcf_default_folder = os.path.join(varank_processing_folder, 'VCF', 'DEFAULT')

		for other_varank_processing_vcf_file in other_varank_processing_vcf_files:
			shutil.move(other_varank_processing_vcf_file, varank_processing_vcf_default_folder)

		varank_processing_vcf_files_default_folder = glob.glob(os.path.join(varank_processing_vcf_default_folder, '*.final.vcf*'))
		for varank_processing_vcf_file_default_folder in varank_processing_vcf_files_default_folder:
			subprocess.run(['rsync', '-rp', varank_processing_vcf_file_default_folder, varank_processing_folder])

	if os.path.isdir(varank_processing_tsv_folder):
		with open(os.path.join(varank_processing_tsv_folder,'TSVDeleting.txt'), mode='a'): pass
		shutil.rmtree(varank_processing_tsv_folder)
		print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') +
			  '] VARANK: Cleaning old TSV folder from ' + varank_processing_folder)
		os.mkdir(varank_processing_tsv_folder)
	else:
		os.mkdir(varank_processing_tsv_folder)


def configfile_family_manager(args, run_platform_application):
	configfile_shelter = os.path.join(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_CONFIG'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_MODULE_NAME'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_SUBMODULE_NAME'], 'configfiles', 'configfile.' + run_platform_application)
	default_configfile_shelter = os.path.join(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_CONFIG'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_MODULE_NAME'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_SUBMODULE_NAME'], 'configfiles', 'configfile.default')
	configfile_shelter_tmp = os.path.join(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_CONFIG'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_MODULE_NAME'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_SUBMODULE_NAME'], 'configfiles', 'configfile.' + run_platform_application + '.tmp')
	default_configfile_shelter_tmp = os.path.join(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_CONFIG'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_MODULE_NAME'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_SUBMODULE_NAME'], 'configfiles', 'configfile.default.tmp')
	varank_family_tracker = os.path.join(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_CONFIG'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_MODULE_NAME'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_SUBMODULE_NAME'], 'configfiles', 'varank_family_tracker.tsv')
	configfiles_folder = os.path.join(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_CONFIG'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_MODULE_NAME'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_SUBMODULE_NAME'], 'configfiles')

	run_platform = run_platform_application.split('.')[0]
	run_application = run_platform_application.split('.')[1]
	families = args.family

	samples_name = []
	samples = glob.glob(os.path.join(args.run, '*', ''))
	for sample in samples:
		sample = os.path.basename(os.path.dirname(sample))
		samples_name.append(sample)

	if os.path.isfile(configfile_shelter):
		fam_conf = configfile_shelter
		tmp_fam_conf = configfile_shelter_tmp
	else:
		print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: There is no specific SampleSheet for this application, adding families to the default SampleSheet')
		fam_conf = default_configfile_shelter
		tmp_fam_conf = default_configfile_shelter_tmp

	all_families = []
	with open(fam_conf, 'r') as read_file:
		for line in read_file.readlines():
			if re.match(r'fam', line):
				line = line.strip()
				samples = re.sub(r'fam[0-9]+: ', '', line)
				samples = samples.split(' ')
				for sample in samples:
					all_families.append(sample)
				number = re.sub('^fam', '', line)
				number = re.sub(':.*', '', number)
			else:
				number = 0

	existing_members_list = []
	member_exist = False

	for family in families:
		family = re.sub(',', ' ', family)
		family = re.sub('\[', '', family)
		family = re.sub('\]', '', family)
		members = family.split(' ')
		if 'fam' in family:
			new_members = []
			for member in members:
				new_member = re.sub(r'^fam[0-9]:', '', member)
				new_members.append(new_member)
			members = new_members

		for member in members:
			if member not in samples_name:
				error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: ' + member + ' doesn\'t exist in this run, please check its spelling, exiting'
				print(error)
				if args.run:
					error_log_writer(error)
					exit()
				exit()

			if member in all_families:
				print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: ' + member + ' was not added to a new family in configfile because he has already its own family, if you want to add a new member to a family you must precise the family number in a specific application, example [fam2:SGT06]')
				member_exist = True
				existing_members_list.append(member)

	if member_exist is True:
		error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] The following samples are already in the configfile and have their own family, please check your' \
				  ' SampleSheet, you can launch VaRank analysis manually with appropriate families : ', \
				existing_members_list
		print(error)
		if args.run:
			error_log_writer(str(error))
			exit()
		exit()

	elif member_exist is False:
		for family in families:
			is_new_member = False
			family = re.sub(',', ' ', family)
			family = re.sub('\[', '', family)
			family = re.sub('\]', '', family)
			members = family.split(' ')
			if 'fam' in family:
				new_members = ''
				for member in members:
					new_member = re.sub(r'^fam[0-9]:', '', member)
					new_member_family = re.sub(r':.*', '', member)
					new_members = new_member + ' ' + new_members
				is_new_member = True

			if is_new_member is False:
				with open(fam_conf, 'a+') as write_file, open(varank_family_tracker, 'a+') as write_file2:
					number = int(number) + 1
					fam_number = 'fam' + str(number)
					fam = fam_number + ': ' + family
					write_file.write(fam + '\n')
					write_file2.write(fam_number + '\t' + family + '\t' + run_platform + '\t' + run_application + '\t' + datetime.now().strftime('%y%y-%m-%d') + '\t' + 'creation' + '\n')
					print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S')
						  + '] VARANK: Added families ' + fam + ' to the ' + fam_conf + ' SampleSheet')

			elif is_new_member is True:
				with open(fam_conf, 'r') as read_file, open(tmp_fam_conf, 'w') as write_file, \
						open(varank_family_tracker, 'a+') as write_file2:
					for line in read_file.readlines():
						line = line.strip()
						if re.match(new_member_family, line):
							line = re.sub(line, line + ' ' + new_members, line)
							family = re.sub(r'fam[0-9]+: ', '', line)
							write_file.write(line + '\n')
							write_file2.write(new_member_family + '\t' + family + '\t' + run_platform + '\t' + run_application + '\t' + datetime.now().strftime('%y%y-%m-%d') + '\t' + 'update' + '\n')
							print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: Family ' + line + ' was updated in ' + fam_conf + ' SampleSheet')
						else:
							write_file.write(line + '\n')

				os.remove(fam_conf)
				os.rename(tmp_fam_conf, fam_conf)
				os.chmod(fam_conf, 0o777)

	configfile_checker(configfiles_folder, args)


def varank_processing_folder_configfile_manager(run_platform_application, varank_processing_folder, args):
	if run_platform_application == 'default':
		configfile_shelter = os.path.join(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_CONFIG'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_MODULE_NAME'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_SUBMODULE_NAME'], 'configfiles', 'configfile.default')
		used_configfile = os.path.join(varank_processing_folder, 'configfile.default')
		renamed_used_configfile = os.path.join(varank_processing_folder, 'configfile')
	else:
		configfile_shelter = os.path.join(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_CONFIG'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_MODULE_NAME'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_SUBMODULE_NAME'], 'configfiles', 'configfile.' + run_platform_application)
		used_configfile = os.path.join(varank_processing_folder, 'configfile.' + run_platform_application)
		configfile_shelter_default = os.path.join(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_CONFIG'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_MODULE_NAME'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_SUBMODULE_NAME'], 'configfiles', 'configfile.default')
		used_default_configfile = os.path.join(varank_processing_folder, 'configfile.default')
		renamed_used_configfile = os.path.join(varank_processing_folder, 'configfile')

	if os.path.isfile(configfile_shelter) and not os.path.isfile(renamed_used_configfile):
		subprocess.run(['rsync', '-rp', configfile_shelter, varank_processing_folder])
		os.rename(used_configfile, renamed_used_configfile)
		print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: Configfile in ' + varank_processing_folder + ' was synced from configfile shelter')
	elif os.path.isfile(configfile_shelter) and os.path.isfile(renamed_used_configfile) and not args.varank_folder:
		subprocess.run(['rsync', '-rp', configfile_shelter, varank_processing_folder])
		os.rename(used_configfile, renamed_used_configfile)
		print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: Configfile in ' + varank_processing_folder + ' was updated from configfile shelter')
	elif not os.path.isfile(configfile_shelter) and os.path.isfile(renamed_used_configfile):
		print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: Using existing configfile in ' + varank_processing_folder)
	elif not os.path.isfile(configfile_shelter) and not os.path.isfile(renamed_used_configfile):
		subprocess.run(['rsync', '-rp', configfile_shelter_default, varank_processing_folder])
		os.rename(used_default_configfile, renamed_used_configfile)
		print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: Using default configfile from configfile shelter')


def varank_launcher(varank_processing_folder):
	logfile = os.path.join(varank_processing_folder, 'VaRank.log')
	varank_bin = os.path.join(os.environ['VARANK'], 'bin', 'VaRank')
	print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: Generating TSV files, you can follow the progress in ' + logfile)
	with open(logfile, 'w') as f:
		subprocess.call([varank_bin, '-vcfdir', varank_processing_folder, '-alamutHumanDB', 'hg19', '-SamVa', '"yes"', '-AlamutProcesses', os.environ['DOCKER_VARANK_USED_THREADS']], stdout=f, stderr=subprocess.STDOUT,universal_newlines=True)
	os.chmod(logfile, 0o777)


def varank_processing_folder_cleaner(varank_processing_tsv_folder, varank_processing_folder, args):
	tsv_files = glob.glob(os.path.join(varank_processing_folder, '*tsv'))
	vcf_files = glob.glob(os.path.join(varank_processing_folder, '*final.vcf.gz'))

	for file in tsv_files:
		os.chmod(file, 0o777)
		shutil.move(file, varank_processing_tsv_folder)

	with open(os.path.join(varank_processing_tsv_folder,'TSVCopyComplete.txt'), mode='a'): pass

	for file in vcf_files:
		os.remove(file)

	generated_tsv_checker(varank_processing_tsv_folder, args)


def generate_non_redundant(varank_processing_tsv_folder, platform_application):
	all_variants_ranking_by_var = os.path.join(varank_processing_tsv_folder, '*_allVariants.rankingByVar.tsv')
	all_variants_ranking_by_var_files = glob.glob(all_variants_ranking_by_var)
	all_variants_ranking_by_var_concat = os.path.join(varank_processing_tsv_folder, 'temp_concat_file.tsv')
	non_redundant_unwanted_columns_default = []
	non_redundant_unwanted_columns = []
	non_redundant_file = os.path.join(varank_processing_tsv_folder, 'Non_Redondant.tsv')

	with open(os.path.join(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_CONFIG'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_MODULE_NAME'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_SUBMODULE_NAME'], 'configfiles', 'non_redundant_config.txt'), 'r') as read_file:
		for line in read_file.readlines():
			line = line.strip()
			if re.match('UNWANTED_VARANK_COLUMN.default', line):
				line = re.sub(r'^.*?=', '', line)
				non_redundant_unwanted_columns_default = line.split(',')
			elif re.match('UNWANTED_VARANK_COLUMN.' + platform_application, line):
				line = re.sub(r'^.*?=', '', line)
				non_redundant_unwanted_columns = line.split(',')
	with open(all_variants_ranking_by_var_concat, 'w') as write_file:
		for all_variants_ranking_by_var_file in all_variants_ranking_by_var_files:
			with open(all_variants_ranking_by_var_file, 'r') as read_file:
				for line in read_file.readlines():
					if not re.match('##', line):
						write_file.write(line)
						
	with open(all_variants_ranking_by_var_concat, 'r'):
		df = pandas.read_csv(all_variants_ranking_by_var_concat, sep='\t', header=[0], low_memory=False, keep_default_na=False)
		if len(non_redundant_unwanted_columns) == 0:
			df = df.drop(non_redundant_unwanted_columns_default, axis=1)
		else:
			df = df.drop(non_redundant_unwanted_columns, axis=1)
		df = df[df.variantID != list(df.columns)[0]]
		if re.match(r'DIAG\.WES_.*', platform_application):
			df['gnomadAltFreq_popmax'] = df['gnomadAltFreq_popmax'].replace({',': '.'}, regex=True)
			index_names = df[df['gnomadAltFreq_popmax'].astype(float) < 0.005].index
			df = df.drop(index_names)
			index_names = df[df['gnomadHomCount_all'].astype(int) < 2].index
			df = df.drop(index_names)
		df = df.drop_duplicates(subset='variantID')
		df.to_csv(non_redundant_file, sep='\t', index=False)

	os.chmod(non_redundant_file, 0o777)
	non_redundant_checker(non_redundant_file)
	os.remove(all_variants_ranking_by_var_concat)
	print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: Non-redondant file was generated')


def run_varank_results_provider(varank_processing_tsv_folder, run_depository, run_repository, varank_processing_folder, run_platform, run_application):
	varank_processing_folder_tsv_files = glob.glob(os.path.join(varank_processing_tsv_folder, '*'))
	run_repository_sample_folder_list = glob.glob(os.path.join(run_repository, '*', ''))

	varank_processing_folder_logfile = os.path.join(varank_processing_folder, 'VaRank.log')
	varank_logfile_depository = os.path.join(run_depository, 'VaRank.log')
	varank_logfile_repository = os.path.join(run_repository, 'VaRank.log')

	varank_renamed_logfile_depository = os.path.join(run_depository, 'VARANK.' + datetime.now().strftime('%y%m%d-%H%M%S') + '.VaRank.log')
	varank_renamed_logfile_repository = os.path.join(run_repository, 'VARANK.' + datetime.now().strftime('%y%m%d-%H%M%S') + '.VaRank.log')

	non_redundant_file_depository = os.path.join(run_depository, 'Non_Redondant.tsv')
	non_redundant_file_repository = os.path.join(run_repository, 'Non_Redondant.tsv')

	non_redundant_renamed_file_depository = os.path.join(run_depository, 'VARANK.' + datetime.now().strftime('%y%m%d-%H%M%S') + '.Non_Redondant.tsv')
	non_redundant_renamed_file_repository = os.path.join(run_repository, 'VARANK.' + datetime.now().strftime('%y%m%d-%H%M%S') + '.Non_Redondant.tsv')

	subprocess.run(['rsync', '-rp', varank_processing_folder_logfile, run_repository])
	subprocess.run(['rsync', '-rp', varank_processing_folder_logfile, run_depository])

	os.rename(varank_logfile_depository, varank_renamed_logfile_depository)
	os.rename(varank_logfile_repository, varank_renamed_logfile_repository)

	for tsv_file in reversed(varank_processing_folder_tsv_files):
		if os.path.basename(tsv_file).startswith('Non_Redondant'):
			subprocess.run(['rsync', '-rp', tsv_file, run_repository])
			subprocess.run(['rsync', '-rp', tsv_file, run_depository])
			varank_processing_folder_tsv_files.remove(tsv_file)
		else:
			for sample_folder in run_repository_sample_folder_list:
				sample_folder = (os.path.basename(os.path.dirname(sample_folder)))
				run_repository_sample_folder_varank_folder = os.path.join(run_repository, os.path.basename(sample_folder), 'VARANK', '')
				run_depository_sample_folder_varank_folder = os.path.join(run_depository, os.path.basename(sample_folder), 'VARANK', '')
				if os.path.basename(tsv_file).startswith('fam'):
					test = re.search(r'(?<=_)(.*?)(?=_[a-zA-Z0-9\.]+.tsv)', os.path.basename(tsv_file)).group(0)
					if test == os.path.basename(sample_folder):
						subprocess.run(['rsync', '-rp', tsv_file, run_repository_sample_folder_varank_folder])
						subprocess.run(['rsync', '-rp', tsv_file, run_depository_sample_folder_varank_folder])
				else:
					subprocess.run(['rsync', '-rp', tsv_file, run_repository_sample_folder_varank_folder])
					subprocess.run(['rsync', '-rp', tsv_file, run_depository_sample_folder_varank_folder])

	os.rename(non_redundant_file_depository, non_redundant_renamed_file_depository)
	os.rename(non_redundant_file_repository, non_redundant_renamed_file_repository)

	varank_complete_file = os.path.join(run_repository, 'VARANKComplete.txt')
	os.rename(os.path.join(run_repository, 'VARANKRunning.txt'), varank_complete_file)
	os.chmod(varank_complete_file, 0o777)

	varank_copy_complete_file = os.path.join(run_depository, 'VARANKCopyComplete.txt')
	subprocess.run(['rsync', '-rp', varank_complete_file, run_depository])
	os.rename(os.path.join(run_depository, 'VARANKComplete.txt'), varank_copy_complete_file)


def launch_run_analysis(args):
	run_depository = run_path_generator(args)[0]
	run_repository = run_path_generator(args)[1]
	if run_repository.endswith('/'):
		run_repository = run_repository[:-1]
	if run_depository.endswith('/'):
		run_depository = run_depository[:-1]

	pattern = pattern_generator(args)
	varank_running_file = os.path.join(run_repository, 'VARANKRunning.txt')
	with open(varank_running_file, 'w') as write_file:
		write_file.write("VaRank is running !")

	path_checker(args, run_depository)
	pattern_checker(pattern, run_repository, args)

	run_application = os.path.basename(os.path.dirname(run_repository))
	run_platform = os.path.basename(os.path.dirname(os.path.dirname(run_repository)))
	run_platform_application = run_platform + '.' + run_application
	varank_processing_folder = os.path.join(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_SERVICES'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_MODULE_NAME'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_SUBMODULE_NAME'], 'Archives', run_platform, run_application)
	varank_processing_tsv_folder = os.path.join(varank_processing_folder, 'TSV')
	varank_processing_run_vcf_folder = os.path.join(varank_processing_folder, 'VCF', os.path.basename(run_repository))
	results = run_repository

	# Steps to produce results
	vcf_synchronizer(pattern, run_repository, varank_processing_run_vcf_folder, varank_processing_folder)
	vcf_synchronization_checker(varank_processing_folder)
	varank_processing_folder_initializer(varank_processing_tsv_folder, varank_processing_folder)
	if args.family:
		configfile_family_manager(args, run_platform_application)
	varank_processing_folder_configfile_manager(run_platform_application, varank_processing_folder, args)
	varank_launcher(varank_processing_folder)
	logfile_checker(varank_processing_folder, args, run_depository, run_repository)
	varank_processing_folder_cleaner(varank_processing_tsv_folder, varank_processing_folder, args)
	generate_non_redundant(varank_processing_tsv_folder, run_platform_application)
	run_varank_results_provider(varank_processing_tsv_folder, run_depository, run_repository, varank_processing_folder, run_platform, run_application)
	print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: Your analysis ended well, you can check results in ' + results)


def launch_folder_analysis(args):
	folder = args.varank_folder
	if folder.endswith('/'):
		folder = folder[:-1]

	if folder.startswith(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_SERVICES']):
		folder_platform = os.path.basename(os.path.dirname(folder))
		folder_application = os.path.basename(folder)
		varank_application_paths = glob.glob(os.path.join(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_SERVICES'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_MODULE_NAME'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_SUBMODULE_NAME'], 'Archives', '*', '*', ''))
		check = True

		for varank_application_path in varank_application_paths:
			varank_application_path = varank_application_path.split('/')
			if not folder_platform == varank_application_path[6] and not folder_application == varank_application_path[7]:
				check = False
			if folder_platform == varank_application_path[6] and folder_application == varank_application_path[7]:
				check = True
				break


		if check is False:
			print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: You are not allowed to launch a VaRank analysis from here ' + folder + ', exiting')
			exit()
		folder_platform_application = folder_platform + '.' + folder_application
		varank_processing_folder = folder
		varank_processing_tsv_folder = os.path.join(varank_processing_folder, 'TSV')
		results = folder

		# Steps to produce results
		varank_processing_folder_initializer(varank_processing_tsv_folder, varank_processing_folder)
		varank_processing_folder_configfile_manager(folder_platform_application, varank_processing_folder)
		varank_launcher(varank_processing_folder)
		logfile_checker(varank_processing_folder, args, 'None', 'None')
		varank_processing_folder_cleaner(varank_processing_tsv_folder, varank_processing_folder, args)
		generate_non_redundant(varank_processing_tsv_folder, folder_platform_application)
		print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: Your analysis ended well, you can check results in ' + results)
	else:
		varank_processing_folder = folder
		varank_processing_tsv_folder = os.path.join(varank_processing_folder, 'TSV')
		folder_platform_application = 'default'
		results = folder

		# Steps to produce results
		varank_processing_folder_initializer(varank_processing_tsv_folder, varank_processing_folder)
		varank_processing_folder_configfile_manager(folder_platform_application, varank_processing_folder, args)
		varank_launcher(varank_processing_folder)
		logfile_checker(varank_processing_folder, args, 'None', 'None')
		varank_processing_folder_cleaner(varank_processing_tsv_folder, varank_processing_folder, args)
		generate_non_redundant(varank_processing_tsv_folder, "default")
		print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: Your analysis ended well, you can check results in ' + results)


def run_path_generator(args):
	if (args.run).endswith('/'):
		run_repository = args.run[:-1]
	else:
		run_repository = args.run
	run_depository = re.sub(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_REPOSITORY'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_DEPOSITORY'], run_repository)
	return run_depository, run_repository


def pattern_generator(args):
	# Defining used pattern
	pattern = []
	check = False

	if args.pattern is None:
		pattern.append('*/STARK/*.reports/*.final.vcf.gz')
		return pattern
	else:
		for element in args.pattern:
			if element == '*/STARK/*.reports/*.final.vcf.gz':
				check = True
				pattern.insert(0, '*/STARK/*.reports/*.final.vcf.gz')
			else:
				pattern.append(element)

		if check is not True:
			pattern.insert(0, '*/STARK/*.reports/*.final.vcf.gz')

		return pattern


def error_log_writer(error):
	args = parse_args()
	run_repository = run_path_generator(args)[1]
	print(run_repository)
	error = error + '\n'
	with open(os.path.join(run_repository, 'VARANKRunning.txt'), 'a+') as write_file:
		write_file.write(error)
	os.chmod(os.path.join(run_repository, 'VARANKRunning.txt'), 0o777)


def generated_tsv_checker(varank_processing_tsv_folder, args):
	tsv_files = glob.glob(os.path.join(varank_processing_tsv_folder, '*tsv'))

	if len(tsv_files) <= 2:
		error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: No TSV files were generated during the analysis, exiting'
		print(error)
		if args.run:
			error_log_writer(error)
			exit()
		exit()
	for tsv_file in tsv_files:
		count_family_barcode = 0
		if os.stat(tsv_file).st_size == 0:
			error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: There are empty TSV files in ' + varank_processing_tsv_folder + ', exiting'
			print(error)
			if args.run:
				error_log_writer(error)
				exit()
			exit()
		if re.match(r'.*rankingByVar.tsv', tsv_file):
			with open(tsv_file, 'r') as read_file:
				for line in read_file.readlines():
					if re.match(r'## FamilyBarcode:.*', line):
						count_family_barcode += 1

		if count_family_barcode > 1:
			error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: Results in *rankingByVar are double for a sample, you need to delete Alamut folder and launch the run again, exiting'
			print(error)
			if args.run:
				error_log_writer(error)
				exit()
			exit()


def vcf_synchronization_checker(varank_processing_folder):
	vcf_folder = os.path.join(varank_processing_folder, 'VCF')

	# Checking if vcf folder exist
	if not os.path.isdir(vcf_folder):
		error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S'), '] VARANK: There is no VCF folder in' + varank_processing_folder + ', exiting'
		print(error)
		error_log_writer(error)
		exit()
	# Checking if there are empty vcf files in vcf folder
	else:
		vcf_files = subprocess.check_output(['find', varank_processing_folder, '-name', '*.final.vcf.gz'])[:-1].decode('utf-8')
		vcf_files_list = vcf_files.split('\n')
		for vcf_file in vcf_files_list:
			if os.stat(vcf_file).st_size == 0:
				error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S'), '] VARANK: There is no VCF files in ' + vcf_folder + ', exiting'
				print(error)
				error_log_writer(error)
				exit()


def path_checker(args, path):
	# Checking if folders exists
	run_repository_folder = os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_REPOSITORY']
	run_depository_folder = os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_DEPOSITORY']

	if not os.path.isdir(path) and path.startswith(run_depository_folder):
		print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: Specified ' + path + ' folder doesn\'t exist, maybe it was sent to Archives ? Creating new folder with all sample subfolders')
		run_repository_sample_folder_list = glob.glob(os.path.join(args.run, '*', ''))
		os.makedirs(path, 0o775)
		for sample_folder in run_repository_sample_folder_list:
			sample_folder = re.sub(run_repository_folder, run_depository_folder, sample_folder)
			if not os.path.isdir(sample_folder):
				os.mkdir(sample_folder, 0o755)

		# Checking if folders are empty
		if len(os.listdir(path)) == 0:
			print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: Specified ' + path + ' folder is empty, exiting')
			exit()


def pattern_checker(pattern, run_repository, args):
	# Checking if there are file(s) respecting the provided pattern, in the case of STARK run analysis
	if args.run:
		check1 = False
		check2 = False
		
	print(run_repository)
	if run_repository.endswith("/"):
		run_repository = run_repository[:-1]
		print(run_repository)

	for element in pattern:
		vcf_files = glob.glob(os.path.join(run_repository, element))
		if len(vcf_files) == 0 and element != '*/*.final.vcf.gz':
			print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: There is no vcf files with the following pattern in this run : ' + element)
		elif len(vcf_files) == 0 and element == '*/*.final.vcf.gz':
			error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: There is no vcf files with the STARK analysis pattern in this run : ' + element
			error_log_writer(error)
			check1 = True

		for vcf_file in vcf_files:
			if os.stat(vcf_file).st_size == 0:
				error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: The following vcf file is empty : ' + vcf_file
				error_log_writer(error)
				check2 = True

	if check1 is True:
		print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: There is no vcf files with the following pattern in this run, exiting')
		if check2 is False:
			exit()
	if check2 is True:
		print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: There are empty vcf files in your run, exiting')
		exit()


def environment_checker(args):
	alamut_license_checker(args)

	# Checking if databases are missing
	if not os.path.isfile(os.path.join(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_DATABASES'], 'OMIMannotations', 'current', 'OMIMannotations.tsv.gz')):
		error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: OMIM database is missing in the ExtAnn VaRank folder, please check docker-compose.yml file, exiting'
		print(error)
		if args.run:
			error_log_writer(error)
			exit()
		exit()
	if not os.path.isfile(os.path.join(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_DATABASES'], 'alamutDB', 'current', 'alamut_db')):
		error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: ALAMUT database is missing in the database folder, please check docker-compose.yml file, exiting'
		print(error)
		if args.run:
			error_log_writer(error)
			exit()
		exit()
	# Checking if databases are empty
	if os.stat(os.path.join(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_DATABASES'], 'OMIMannotations', 'current', 'OMIMannotations.tsv.gz')).st_size == 0:
		error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: OMIM database is empty in the ExtAnn VaRank folder, please check docker-compose.yml file, exiting'
		print(error)
		if args.run:
			error_log_writer(error)
			exit()
		exit()
	if os.stat(os.path.join(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_DATABASES'], 'alamutDB', 'current', 'alamut_db')).st_size == 0:
		error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: ALAMUT database is empty in the database folder, please check docker-compose.yml file, exiting'
		print(error)
		if args.run:
			error_log_writer(error)
			exit()
		exit()

	# Checking if mounted volumes are missing
	if not os.path.isdir(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_REPOSITORY']):
		error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: STARK repository folder is missing inside the container, please check your docker configuration, exiting'
		print(error)
		if args.run:
			error_log_writer(error)
			exit()
		exit()
	if not os.path.isdir(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_DEPOSITORY']):
		error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: STARK depository folder is missing inside the container, please check your docker configuration, exiting'
		print(error)
		if args.run:
			error_log_writer(error)
			exit()
		exit()
	if not os.path.isdir(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_CONFIG']):
		error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: VARANK config folder is missing inside the container, please check your docker configuration, exiting'
		print(error)
		if args.run:
			error_log_writer(error)
			exit()
		exit()
	# Checking if mounted volumes are empty
	if len(os.listdir(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_REPOSITORY'])) == 0:
		error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: STARK repository folder is empty inside the container, please check your docker configuration, exiting'
		print(error)
		if args.run:
			error_log_writer(error)
			exit()
		exit()
	if len(os.listdir(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_CONFIG'])) == 0:
		error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: VARANK config folder is empty inside the container, please check your docker configuration, exiting'
		print(error)
		if args.run:
			error_log_writer(error)
			exit()
		exit()


def varank_config_folder_checker(args):
	configfiles_folder = os.path.join(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_CONFIG'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_MODULE_NAME'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_SUBMODULE_NAME'], 'configfiles')
	varank_config = os.path.join(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_CONFIG'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_MODULE_NAME'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_SUBMODULE_NAME'])
	extanns = os.path.join(configfiles_folder, 'extanns')
	non_redundant_config = os.path.join(configfiles_folder, 'non_redundant_config.txt')
	deleted_samples_config = os.path.join(configfiles_folder, 'deleted_samples_config.tsv')
	extann_config_file = os.path.join(configfiles_folder, 'extann_config_file.tsv')
	varank_family_tracker = os.path.join(configfiles_folder, 'varank_family_tracker.tsv')
	configfile_default = os.path.join(configfiles_folder, 'configfile.default')
	varank_config_initialized = os.path.join(varank_config, 'initialized.txt')

	alamut_license = os.path.join(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_CONFIG'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_MODULE_NAME'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_SUBMODULE_NAME'], 'alamut-batch-license', 'alamut-batch.ini')
	alamut_license_folder = os.path.join(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_CONFIG'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_MODULE_NAME'], os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_SUBMODULE_NAME'], 'alamut-batch-license')

	if not os.path.isdir(varank_config):
		os.mkdir(varank_config, 0o777)
	elif os.path.isdir(varank_config):
		os.chmod(varank_config, 0o777)
	if args.configure:
		with open(varank_config_initialized, 'a+') as write_file:
			write_file.write(datetime.now().strftime('%y%y%m%d-%H%M%S') + '\n')
			write_file.write('Do not deleted this file !' + '\n')
			pass
		os.chmod(varank_config_initialized, 0o755)
	if not os.path.isfile(varank_config_initialized):
		error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: For a first launch, you need to launch manually VaRank with --configure option, exiting'
		print(error)
		if args.run:
			error_log_writer(error)
			exit()
		exit()
	if not os.path.isdir(configfiles_folder):
		os.mkdir(configfiles_folder, 0o777)
	elif os.path.isdir(configfiles_folder):
		os.chmod(configfiles_folder, 0o777)
		configfile = glob.glob(os.path.join(configfiles_folder, 'configfile.*'))
		for file in configfile:
			os.chmod(file, 0o777)
	if not os.path.isfile(configfile_default):
		shutil.copy(os.path.join(os.environ['VARANK'], 'configfile'), os.path.join(configfiles_folder, 'configfile.default'))
	elif os.path.isfile(configfile_default):
		os.chmod(configfile_default, 0o777)
	if not os.path.isdir(extanns):
		os.mkdir(extanns, 0o777)
	elif os.path.isdir(extanns):
		os.chmod(extanns, 0o777)
	if not os.path.isfile(non_redundant_config):
		shutil.copy(os.path.join(os.environ['VARANK'], 'config', 'non_redundant_config.txt'), configfiles_folder)
		os.chmod(non_redundant_config, 0o777)
	elif os.path.isfile(non_redundant_config):
		os.chmod(non_redundant_config, 0o777)
	if not os.path.isfile(deleted_samples_config):
		shutil.copy(os.path.join(os.environ['VARANK'], 'config', 'deleted_samples_config.tsv'), configfiles_folder)
		os.chmod(deleted_samples_config, 0o777)
	elif os.path.isfile(deleted_samples_config):
		os.chmod(deleted_samples_config, 0o777)
	if not os.path.isfile(extann_config_file):
		shutil.copy(os.path.join(os.environ['VARANK'], 'config', 'extann_config_file.tsv'), configfiles_folder)
		os.chmod(extann_config_file, 0o777)
	elif os.path.isfile(extann_config_file):
		os.chmod(extann_config_file, 0o777)
	if not os.path.isfile(varank_family_tracker):
		shutil.copy(os.path.join(os.environ['VARANK'], 'config', 'varank_family_tracker.tsv'), configfiles_folder)
		os.chmod(varank_family_tracker, 0o777)
	elif os.path.isfile(varank_family_tracker):
		os.chmod(varank_family_tracker, 0o777)

	if not os.path.isdir(alamut_license_folder):
		os.mkdir(alamut_license_folder, 0o777)
		if not os.path.isfile(alamut_license):
			error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: Using default alamut license file within the installation folder, you need to add your specific license in ' + alamut_license + ', exiting'
			subprocess.run(['rsync', '-rp', os.path.join(os.environ['ALAMUT'], 'alamut-batch.ini'), alamut_license])
			print(error)
			if args.run:
				error_log_writer(error)
				exit()
			exit()
	elif os.path.isfile(alamut_license):
		subprocess.run(['rsync', '-rp', alamut_license, os.path.join(os.environ['ALAMUT'], 'alamut-batch.ini')])
		os.chmod(alamut_license, 0o777)
	configfile_checker(configfiles_folder, args)


def configfile_checker(configfiles_folder, args):
	configfiles_list = glob.glob(os.path.join(configfiles_folder, 'configfile*'))
	extann_config_file = os.path.join(configfiles_folder, 'extann_config_file.tsv')
	extann_directory = os.path.join(configfiles_folder, 'extanns')
	vcf_fields = False

	for configfile in configfiles_list:
		family_list = []
		tmp_configfile = os.path.join(os.path.dirname(configfile), 'tmp' + os.path.basename(configfile))
		with open(configfile, 'r') as read_file:
			for line in read_file.readlines():
				line = line.strip()
				if re.match(r'^-vcfFields:', line):
					vcf_fields = True

		with open(tmp_configfile, 'w') as write_file, open(configfile, 'r') as read_file, \
				open(extann_config_file, 'r') as read_file2:
			for line in read_file.readlines():
				writted = False
				line = line.strip()

				if re.match(r'-vcfInfo:|#-vcfInfo:', line):
					if re.match(r'#-vcfInfo:', line):
						line = re.sub('#', '', line)
					if re.search('no', line):
						line = re.sub('no', 'yes', line)
					if vcf_fields is False:
						write_file.write(line + '\n' + '-vcfFields: "FindByPipelines GenotypeConcordance POOL_F_Depth POOL_M_Depth POOL_F_base_counts POOL_M_base_counts BARCODE trio_variant_type"' + '\n')
						writted = True
					elif vcf_fields is True:
						write_file.write(line + '\n')
						writted = True
					vcf_fields = True

				if re.match(r'-metrics:|#-metrics:', line):
					if re.match(r'#-metrics:', line):
						line = re.sub('#', '', line)
					if re.search('us', line):
						line = re.sub('us', 'fr', line)
					write_file.write(line + '\n')
					writted = True

				if re.match(r'#-uniprot:', line):
					line = re.sub('#', '', line)
					write_file.write(line + '\n')
					writted = True

				if re.match(r'#-refseq:', line):
					line = re.sub('#', '', line)
					write_file.write(line + '\n')
					writted = True

				if re.match(r'fam[0-9]+:', line):
					fam = line.split(' ')
					fam = fam[1:]
					for id in fam:
						family_list.append(id)

				if re.match(r'-extann:|#-extann:', line):
					extann_files = line.split('"')
					for extann_file in reversed(extann_files):
						if not extann_file.startswith('/'):
							extann_files.remove(extann_file)

					if os.path.basename(configfile) != 'configfile.default':
						if len(extann_files) != 0:
							for extann_file in extann_files:
								extann_file_list = extann_file.split(' ')
						else:
							extann_file_list = []

						for extann_file in reversed(extann_file_list):
							if re.search(r'^.*$', extann_file):
								extann_file_list.remove(extann_file)

						next(read_file2, 1)
						for line2 in read_file2.readlines():
							line2 = line2.strip()
							line2 = line2.split('\t')
							extann_platform_application = line2[0] + '.' + line2[1]
							extann_name = os.path.join(extann_directory, line2[2])
							extann_used = line2[3]
							configfile_platform_application = os.path.basename(configfile).split('.')[1] + '.' + os.path.basename(configfile).split('.')[2]

							if extann_platform_application == configfile_platform_application and extann_used == 'yes':
								extann_file_list.append(extann_name)

						if len(extann_file_list) == 0 and args.configure:
							print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: Not using any additional extann(s) file(s) for ' + configfile)
						for extann_file in reversed(extann_file_list):
							if os.path.isfile(extann_file) and len(extann_file_list) != 0:
								continue
							elif not os.path.isfile(extann_file) and len(extann_file_list) != 0:
								extann_file_list.remove(extann_file)
								error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: Not found additional extann(s) file(s) at ' + extann_file + 'directory, please check your configfile, exiting'
								print(error)
								if args.run:
									error_log_writer(error)
									exit()
								exit()

						if len(extann_file_list) == 0:
							line = '#-extann:\t\t""'
						else:
							line = '-extann:\t\t"' + " ".join(extann_file_list) + '"'
						write_file.write(line + '\n')
						writted = True

				if re.match(r'#-proxyUser:', line):
					if re.search(r'"[a-zA-Z0-9]+"', line):
						line = re.sub('"[a-zA-Z0-9]+"', '"' + os.environ['PROXY_USER'] + '"', line)
					else:
						line = re.sub('""', '"' + os.environ['PROXY_USER'] + '"', line)
					line = re.sub('#', '', line)
					write_file.write(line + '\n')
					writted = True
				
				if re.match(r'#-proxyPasswd:', line):
					if re.search(r'"[a-zA-Z0-9]+"', line):
						line = re.sub('"[a-zA-Z0-9]+"', '"' + os.environ['PROXY_PASSWD'] + '"', line)
					else:
						line = re.sub('""', '"' + os.environ['PROXY_PASSWD'] + '"', line)
					line = re.sub('#', '', line)
					write_file.write(line + '\n')
					writted = True
				
				if re.match(r'#-proxyServer:', line):
					if re.search(r'"[a-zA-Z0-9]+"', line):
						line = re.sub('"[a-zA-Z0-9]+"', '"' + os.environ['PROXY_SERVER'] + '"', line)
					else:
						line = re.sub('""', '"' + os.environ['PROXY_SERVER'] + '"', line)
					line = re.sub('#', '', line)
					write_file.write(line + '\n')
					writted = True
				
				if re.match(r'#-proxyPort:', line):
					if re.search(r'"[a-zA-Z0-9]+"', line):
						line = re.sub('"[a-zA-Z0-9]+"', '"' + os.environ['PROXY_PORT'] + '"', line)
					else:
						line = re.sub('""', '"' + os.environ['PROXY_PORT'] + '"', line)
					line = re.sub('#', '', line)
					write_file.write(line + '\n')
					writted = True
					
				if re.match(r'-proxyUser:', line):
					if re.search(r'"[a-zA-Z0-9]+"', line):
						line = re.sub('"[a-zA-Z0-9]+"', '"' + os.environ['PROXY_USER'] + '"', line)
					else:
						line = re.sub('""', '"' + os.environ['PROXY_USER'] + '"', line)
					write_file.write(line + '\n')
					writted = True
				
				if re.match(r'-proxyPasswd:', line):
					if re.search(r'"[a-zA-Z0-9]+"', line):
						line = re.sub('"[a-zA-Z0-9]+"', '"' + os.environ['PROXY_PASSWD'] + '"', line)
					else:
						line = re.sub('""', '"' + os.environ['PROXY_PASSWD'] + '"', line)
					write_file.write(line + '\n')
					writted = True
				
				if re.match(r'-proxyServer:', line):
					if re.search(r'"[a-zA-Z0-9]+"', line):
						line = re.sub('"[a-zA-Z0-9]+"', '"' + os.environ['PROXY_SERVER'] + '"', line)
					else:
						line = re.sub('""', '"' + os.environ['PROXY_SERVER'] + '"', line)
					write_file.write(line + '\n')
					writted = True
				
				if re.match(r'-proxyPort:', line):
					if re.search(r'"[a-zA-Z0-9]+"', line):
						line = re.sub('"[a-zA-Z0-9]+"', '"' + os.environ['PROXY_PORT'] + '"', line)
					else:
						line = re.sub('""', '"' + os.environ['PROXY_PORT'] + '"', line)
					write_file.write(line + '\n')
					writted = True

				if writted is False:
					write_file.write(line + '\n')

			# if any(family_list.count(element) > 1 for element in family_list):
			#	  error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + \
			#			  '] VARANK: There are multiple time same patient(s) in your ' + configfile + \
			#			  ', please keep only one, exiting'
			#	  print(error)
			#	  if args.run:
			#		  error_log_writer(error)
			#		  exit()
			#	  exit()

		os.remove(configfile)
		os.rename(tmp_configfile, configfile)
		os.chmod(configfile, 0o777)


def logfile_checker(varank_processing_folder, args, run_depository, run_repository):
	varank_error = True
	logfile = os.path.join(varank_processing_folder, 'VaRank.log')

	with open(logfile, 'r') as read_file:
		for line in read_file.readlines():
			if re.match(r'^\.\.\.VaRank is done with the analysis', line):
				varank_error = False

	if varank_error is True:
		error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: There is an unexpected error in the VaRank analysis, you must check its log, exiting'
		print(error)
		if args.run:
			error_log_writer(error)
			subprocess.run(['rsync', '-rp', logfile, run_depository])
			subprocess.run(['rsync', '-rp', logfile, run_repository])
			exit()
		exit()


def non_redundant_checker(non_redundant_file):
	with open(non_redundant_file, 'r') as read_file:
		df = pandas.read_csv(read_file, sep='\t', header=[0], low_memory=False)
		df.duplicated(subset=['variantID'])
		for check in df.duplicated(subset=['variantID']):
			if check is True:
				error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: There are duplicated variantID in your non-redundant file, exiting'
				print(error)
				error_log_writer(error)
				exit()
				break


def alamut_license_checker(args):
	with open(os.path.join(os.environ['ALAMUT'], 'alamut-batch.ini'), 'r') as read_file:
		for line in read_file:
			line = line.strip()
			if line == 'licenseKey=...':
				error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: There is no licenseKey defined in the alamut-batch.ini file, exiting'
				print(error)
				if args.run:
					error_log_writer(error)
					exit()
				exit()

			if re.search(r'(Institution=[A-Za-z]+)', line):
				print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: licenseKey is valid for user', line.strip('Institution='))
			if re.search('File=', line):
				line = line.strip('File=')
				if line == os.path.join(os.environ['DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_DATABASES'], 'alamutDB', 'current', 'alamut_db'):
					print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: AlamutDB path is valid in the license file')
				else:
					error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: Wrong path for the AlamutDB in the alamut-batch.ini file, exiting'
					print(error)
					if args.run:
						error_log_writer(error)
						exit()
					exit()


def args_checker(args):
	# Checking arguments conflicts
	if args.varank_folder and args.pattern:
		print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: You cannot use --varank_folder and --pattern together, please use only one of them, exiting')
		exit()
	if args.varank_folder and args.allsync18:
		print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: You cannot use --varank_folder and --allsync18 together, please use only one of them, exiting')
		exit()
	if args.varank_folder and args.configure:
		print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: You cannot use --varank_folder and --configure together, please use only one of them, exiting')
		exit()
	if args.pattern and args.configure:
		print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: You cannot use --pattern and --configure together, please use only one of them, exiting')
		exit()
	if args.run and args.configure:
		error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: You cannot use --run and --configure together, please use only one of them, exiting'
		print(error)
		error_log_writer(error)
		exit()
	if args.run and args.allsync18:
		error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: You cannot use --run and --allsync18 together, please use only one of them, exiting'
		print(error)
		error_log_writer(error)
		exit()
	if args.pattern and not args.run and not args.allsync18:
		error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: You cannot use --pattern alone, exiting'
		print(error)
		error_log_writer(error)
		exit()
	if args.varank_folder and args.pattern and args.allsync18:
		print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: You cannot use --varank_folder, --pattern together and --allsync19 together, please use only' ' one of them, exiting')
		exit()
	if args.run and args.allsync18 and args.pattern:
		error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: You cannot use --run, --allsync18 and --pattern together, please use only one of them, exiting'
		print(error)
		error_log_writer(error)
		exit()

	arg_check = False
	args_dict = vars(args).values()
	for value in args_dict:
		if value is not None and value is not False:
			arg_check = True
			break

	if arg_check is False:
		error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: You need to enter a parameter, type -h to see help, exiting'
		print(error)
		if args.run:
			error_log_writer(error)
			exit()
		exit()

	# Checking if specified path are absolute
	if args.varank_folder and not os.path.isabs(args.varank_folder):
		print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: You have not specified an absolute path for --varank_folder argument, exiting')
		exit()

	if args.run and not os.path.isabs(args.run):
		print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: You have not specified an absolute path for --run argument, exiting')
		exit()

	# Checking if specified directories exists
	if args.varank_folder and not os.path.isdir(args.varank_folder):
		print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: Specified VaRank analysis folder doesn\'t exist, exiting')
		exit()
	if args.run and not os.path.isdir(args.run):
		print('[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: Specified VaRank analysis run folder doesn\'t exist, exiting')
		exit()

	# Checking family format
	if args.family:
		for family in args.family:
			if not re.match(r'\[[:a-zA-Z0-9,]+\]', family):
				error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: Error in families format, it must be [sample1,sample2,sample3] [sample4, sample5] for example, exiting'
				print(error)
				if args.run:
					error_log_writer(error)
					exit()
				exit()
			if 'fam' in family:
				family = re.sub('\[', '', family)
				family = re.sub(']', '', family)
				family = re.sub(':[a-zA-Z0-9_]+', '', family)
				family = family.split(',')
				result = False
				if len(family) > 1:
					result = all(elem == family[0] for elem in family)
					if result is False:
						error = '[' + datetime.now().strftime('%y%y%m%d-%H%M%S') + '] VARANK: Error in families format, you are adding samples in pre-existing families you are not allowed to specify many families in one group, you need to declare a new group for each family update'
						print(error)
						if args.run:
							error_log_writer(error)
							exit()
						exit()


def parse_args():
	parser = argparse.ArgumentParser(description='Launch VaRank analysis')
	group = parser.add_mutually_exclusive_group(required=False)
	group.add_argument('-f', '--varank_folder', type=str, help='absolute path to the folder where you want to launch VaRank analysis')
	group.add_argument('-r', '--run', type=str, help='path to STARK run for which you want to launch VaRank analysis')
	parser.add_argument('-p', '--pattern', type=str, nargs='+', help='pattern describing which vcf files to synchronize in STARK folders, actual patterns : */*.final.vcf.gz, */POOL/*.final.vcf.gz. You can use two patterns with space as separator')
	parser.add_argument('-as18', '--allsync18', help='Used to rsync all VCF files from all stark18 project to your VaRank processing directory, no argument required must be used alone', action='store_true')
	parser.add_argument('-fam', '--family', type=str, nargs='+', help='You can pass families to add them to the configfile, the format is the following : -fam [SGT01,SGT02,SGT03] [SGT04,SGT05]')
	parser.add_argument('-cfg', '--configure', help='Use this argument the very first time you launch VaRank to configure all the files and folders', action='store_true')

	return parser.parse_args()


if __name__ == '__main__':
	main()
