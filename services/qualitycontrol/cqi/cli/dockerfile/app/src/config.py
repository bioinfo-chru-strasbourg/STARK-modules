import json
import re 
import os
import sys

#TODO ALREADY add in the listener
#need to export variables ex: export CQI
file = "/home1/BAS/lamouchj/CQI/bin/CQI.json"
#json_raw= raw.readlines()
CQI = os.environ["CQI_SAMPLE"]


def parseJson(file, CQI):
    with open (file) as f:
        jload = json.load(f)
        for item in jload['CQI']:
            #print(item['name'])
            if item['name'] == CQI:
                #print(item['VCF'])
                print("[ #INFO ] CQI " + item['name'] + " detected !")
                CQI_VCF = item['VCF']
                os.putenv("CQI_VCF", "item['VCF']")
                print(CQI_VCF)
                sys.exit()
            else:
                continue
        print(" [ ERROR ] Wrong CQI provided! got "+"{"+CQI+ "}")
        sys.exit()
        return CQI_VCF

parseJson(file, CQI)