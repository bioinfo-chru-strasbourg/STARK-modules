#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function

import os
import re
import subprocess
import glob

from datetime import datetime


def create_container_file(containersFile, run, containerName):
	with open(os.path.join(containersFile, containerName + '.log'), 'w+') as write_file:
		write_file.write('RUN: ' + os.path.basename(run) + '\n')
		write_file.write('FOLDER: ' + run + '\n')
		write_file.write('EXEC_DATE: ' + datetime.now().strftime('%d%m%Y-%H%M%S') + '\n')
		write_file.write('ID: ' + containerName + '\n')


def create_running_file(run, serviceName):
	with open(os.path.join(run, serviceName + 'Running.txt'), 'w+') as write_file:
		write_file.write('# [' + datetime.now().strftime('%d/%m/%Y %H:%M:%S') + '] ' + os.path.basename(run) + ' running with ' + serviceName + '\n')
	os.chmod(os.path.join(run, serviceName + 'Running.txt'), 0o777)


def find_family_tag(run):
	family_list = []
	family_dict = {}
	member_dict = {}
	flipped_family_dict = {}
	sample_folder_list = glob.glob(os.path.join(run, '*', ''))

	for sample_folder in sample_folder_list:
		sample_folder = sample_folder[:-1]
		if sample_folder.endswith('logs') or sample_folder.endswith('CANOES'):
			break
		tag_file = os.path.join(sample_folder, os.path.basename(sample_folder) + '.tag')
		with open(tag_file, 'r') as read_file:
			for tags in read_file:
				tags = tags.strip()
				tags = tags.split('!')
				for tag in tags:
					if tag.startswith('FAMILY#'):
						family_tag = tag.split('#')
						family_dict[os.path.basename(sample_folder)] = family_tag[1]
					if tag.startswith('MEMBER#'):
						family_member = tag.split('#')
						member_dict[os.path.basename(sample_folder)] = family_member[1]

	for family_key in family_dict:
		family_value = family_dict[family_key]
		if family_value not in flipped_family_dict:
			flipped_family_dict[family_value] = [family_key]
		else:
			flipped_family_dict[family_value].append(family_key)

	for flipped_family_key in flipped_family_dict:
		family = '['
		for member in flipped_family_dict[flipped_family_key]:
			family = family + member + ','
		family = family[:-1] + ']'
		family_list.append(family)

	return ' '.join(family_list)


def launch(run, serviceName, containersFile, montage, image, launchCommand, configFile, repository):
	original_umask = os.umask(0o000)
	create_running_file(run, serviceName)
	container_name = serviceName + '_' + os.path.basename(os.path.dirname(run)) + '_' + os.path.basename(run)
	family_tag = find_family_tag(run)
	create_container_file(containersFile, run, container_name)
	if family_tag == '':
		cmd = 'docker run --rm --name=' + container_name + ' ' + montage + ' ' + image + ' -r ' + run
		print(cmd)
		subprocess.call(cmd, shell=True)
	else:
		cmd = 'docker run --rm --name=' + container_name + ' ' + montage + ' ' + image + ' -r ' + run + ' -fam ' + family_tag
		print(cmd)
		subprocess.call(cmd, shell=True)
	os.umask(original_umask)