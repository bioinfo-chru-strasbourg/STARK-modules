#!/usr/bin/python
# -*- coding: latin1 -*-
#

# mispred_ht.py <substs> <seq_file> [<uniprot_acc>]
# substs and optionally uniprot_acc are passed through arguments
# sequences are stored in a tmp file which filename is passed as an argument
# <substs> = <subst>[|<subst>]*
# <subst> = [A-Z][0-9]+[A-Z]

import csv
import os
import re
import signal
import subprocess
import sys
import time

from os.path import exists
from tempfile import NamedTemporaryFile

# debug
debug = False
siftDir = "sift4.0.3b"
tmpdir = "tmp/"

def main():
    save = []
    try:
        if len(sys.argv) < 3:
            print "usage: mispred_ht.py <substs> <seq_file> [<uniprot_acc>]"
            sys.exit(0)

        if (len(sys.argv) == 3):
            uniprot_acc = ''
            seq_file = sys.argv[-1]
            substitutions = sys.argv[-2]
        else:
            uniprot_acc = sys.argv[-1]
            seq_file = sys.argv[-2]
            substitutions = sys.argv[-3]

        fi = open(seq_file)
        aligned_seqs = fi.read()
        fi.close()

        scorerHT = ScorerHT(aligned_seqs, uniprot_acc, substitutions)
        save.append(scorerHT.tmp)
        scorerHT.startxml()
#        scorerHT.agvgd()
#        scorerHT.mapp()
#        scorerHT.pph()
        scorerHT.sift()
        scorerHT.endxml()
    except Exception, errors:
        import time
        import traceback
        errtime = '--- '+ time.ctime(time.time()) +' ---\n'
        errlog = open('error.log', 'a')
        errlog.write(errtime)
        errlog.write(str(save))
        traceback.print_exc(None, errlog)
        errlog.close()

class ScorerHT:
    def __init__(self, aligned_seqs, uniprot_acc, substs):
        self.aligned_seqs = aligned_seqs
        self.uniprot_acc = uniprot_acc

        # tmp filenames
        self.pid = os.getpid()
        prefix = 'ht_%s' % self.pid
        self.tmp = NamedTemporaryFile(prefix=prefix, dir=tmpdir)

        # shared aligned seqs file
        self.aligned_seqs_file = '%s.fasta' % self.tmp.name
        fi = open(self.aligned_seqs_file, 'w')
        fi.write(aligned_seqs)
        fi.close()

        self.substs = {}
        self.positions = {}
        self.maxpos = 0
        self.rejected = []

        # input file
        self.input_file = '%s_input.txt' % self.tmp.name
        fi = open(self.input_file, 'w')
        for string in re.split('\|', substs):
            m = re.match('([A-Z])(\d+)([A-Z])', string)
            if m:
                subst = "%s%s%s" % m.groups()
                self.substs[subst] = m.groups() # {AAnAA => (aa_from, pos, aa_to)}
                fi.write('%s\n' % subst)
                pos = m.group(2)
                self.positions[pos] = 1
                i = int(pos)
                if i > self.maxpos:
                    self.maxpos = i
            elif string:
                self.rejected.append(string)
        fi.close()

        # xml tags
        self.xml_header = '<score source="%s" version="%s">'
        self.xml_status = '<status value="%s" />'

    def pph(self):
        # prog
        self.prog = 'pph'
        self.version = '2.0.22r312'
        self.xmls = []

        flag = True
        if self.uniprot_acc:
            # output file
            output_file = '/home/andre/software/pph2_whpss/%s.pph2.txt' % self.uniprot_acc
            if not exists(output_file):
                self.status = 'no %s in pph2' % self.uniprot_acc
                flag = False
        else:
            self.status = 'no %s in alaref' % uniprot_acc
            flag = False

        if flag:
            # xml tag
            xml_result = '<result aa_from="%s" pos="%s" aa_to="%s" prediction="%s" based_on="%s" effect="%s" pph2_class="%s" pph2_prob="%s" site="%s" region="%s" />'

            # run pph
            self.status = 'ok'
            fi = open(output_file)
            fi.next()
            try:
                for line in fi:
                    row = re.split('[ ]*\t[ ]*', line.rstrip('\n'))
                    pos = row[7] # pos
                    if int(pos) > self.maxpos:
                        break
                    if not pos in self.positions:
                        continue
                    aa_from = row[8] # aa1
                    aa_to = row[9] # aa2
                    subst = "%s%s%s" % (aa_from, pos, aa_to)
                    if not subst in self.substs:
                        continue
                    prediction = row[12] # prediction
                    based_on = row[13] # based_on
                    effect = row[14] # effect
                    pph2_class = row[15] # pph2_class
                    pph2_prob = row[16] # pph2_prob
                    site = row[20] # site
                    region = row[21] # region
                    if prediction == 'unknown':
                        pph2_prob = '0.0'
                    else:
                        float(pph2_prob) # check data
                    self.xmls.append(xml_result % (aa_from, pos, aa_to, prediction, based_on, effect, pph2_class, pph2_prob, site, region))
            except HTException:
                self.status = 'format'

        # xml output
        self.printxml()

    def sift(self):
        # prog
        self.prog = 'sift'
        self.version = '4.0.3'

        # output file
        output_file = '%s_sift_output.txt' % self.tmp.name

        # main cmd
        self.dir = siftDir
        cmd = 'env BLIMPS_DIR=%s/blimps/ %s/bin/info_on_seqs %s %s %s'
        self.cmdline = cmd % (self.dir, self.dir, self.aligned_seqs_file, self.input_file, output_file)
        print self.cmdline
        self.timeoutseconds = 300
        self.xmls = []
        # run

        flag = self.run()
        self.status = 'ok'
        if flag:
            # xml tag
            xml_result = '<result aa_from="%s" pos="%s" aa_to="%s" prediction="%s" weight="%s" median="%s" />'
            fi = open(output_file)
            try:
                for line in fi:
                    row = re.split('\t', line.rstrip('\n'))
                    variant_id = row[0]
                    subst = self.substs.get(variant_id, None)
                    if not subst:
                        continue
                    (aa_from, pos, aa_to) = subst
                    prediction = row[1]
                    if len(row) == 6:
                        weight = row[2]
                        median = row[3]
                        float(weight) # check data
                    else:
                        weight = '0.0'
                        median = '0.0'
                    self.xmls.append(xml_result % (aa_from, pos, aa_to, prediction, weight, median))
            except HTException:
                self.status = 'format'
            fi.close()
        else:
            self.status = 'time out'

        
        os.remove(self.aligned_seqs_file)
        os.remove(self.input_file)
        os.remove(output_file)

        # xml output
        self.printxml()

    def run(self):
        signal.signal(signal.SIGALRM, alarm_handler)
        signal.alarm(self.timeoutseconds)
        try:
            proc = subprocess.Popen(self.cmdline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize =-1, shell=True)
            self.stdout, self.stderr = proc.communicate()
            signal.alarm(0)
            return True
        except HTException:
            return False

    def startxml(self):
        print '<?xml version="1.0"?>'
        if self.rejected:
            print '<scores rejected="%s">' % (','.join(self.rejected))
        else:
            print '<scores>'

    def printxml(self):
        print self.xml_header % (self.prog, self.version)
        print self.xml_status % self.status
        for xml in self.xmls:
            print xml
        if debug:
            print '<cdmline>%s</cdmline>' % self.cmdline
            print '<stdout>%s</stdout>' % self.stdout
            print '<stderr>%s</stderr>' % self.stderr
        print '</score>'

    def endxml(self):
        print '</scores>'

class HTException(Exception):
    pass

def alarm_handler(signum, frame):
    raise HTException

if (__name__ == "__main__"):
    main()
