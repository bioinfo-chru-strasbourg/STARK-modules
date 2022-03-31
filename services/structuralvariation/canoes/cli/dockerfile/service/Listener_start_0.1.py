def find_any_samplesheet(runDir, fromResDir = False):
	"""
	Adapted from runmetrics.py
	
	1) look up recursively all files named SampleSheet.csv in the runDir
	2) check if file path follows an expected samplesheet name and location
		(the latter depends on if we're in a STARK result or repository dir,
		defined by the bool fromResDir)
	3) first correct file path is returned
	"""
	p = subprocess.Popen("find -L "+runDir+" -maxdepth 3 -name *SampleSheet.csv", stdout=subprocess.PIPE, shell=True)
	out = p.stdout.readlines()
	for ss in out:
		ss = ss.decode("utf-8").strip()
		if fromResDir:
			r = re.match(runDir.rstrip("/")+"/(.*)/(.*).SampleSheet.csv", ss)
		else:
			r = re.match(runDir.rstrip("/")+"/(.*)/STARK/(.*).SampleSheet.csv", ss)
		if r is None:
			continue
		elif r.group(1) == r.group(2): #checks if (.*) == (.*)
			return ss
	return "NO_SAMPLESHEET_FOUND"

def createContainerFile(config, run, serviceName):
	if os.path.exists(osj(os.path.dirname(config),"ContainerList"+serviceName+".txt")):
		file = open(osj(os.path.dirname(config),"ContainerList"+serviceName+".txt"), "a")
		file.write(run+"\n")
		file.close()
	else:
		file = open(osj(os.path.dirname(config),"ContainerList"+serviceName+".txt"), "w+")
		file.write(run+"\n")
		file.close()

def createRunningFile(run, serviceName):
	file = open(osj(run,serviceName+"Running.txt"), "w+")
	file.write("# ["+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+"] "+os.path.basename(run)+" running with "+serviceName+"\n")
	file.close()

def getMd5(run):
	runMd5 = hashlib.md5()
	runMd5.update(hashlib.md5(run).hexdigest())
	return runMd5.hexdigest()

def findAnyBed(run):
	p = subprocess.Popen("find -L "+run+" -maxdepth 3 -name *.bed", stdout=subprocess.PIPE, shell=True)
	out = p.stdout.readlines()
	for bed in out:
		bed = bed.decode("utf-8").strip()
		r = re.match(run.rstrip("/")+"/(.*)/STARK/(.*).bed", bed)
		if r is None:
			continue
		elif r.group(1) == r.group(2): #checks if (.*) == (.*)
			return bed
	return "NO_BED_FOUND"

def startService(run, serviceName, serviceDockerImage, annotsvServer, annotsvContainer, genome, config):
	samplesheet = find_any_samplesheet(run)
	assert samplesheet != "NO_SAMPLESHEET_FOUND",\
		"[ERROR] find_any_samplesheet() couldn't find any samplesheet in run"+run+"."
	bed = findAnyBed(run)
	assert bed != "NO_BED_FOUND",\
		"[ERROR] findAnyBed() couldn't find any bed in run"+run+"."
	md5 = getMd5(run)
	containerName = "service-"+serviceName+"-"+md5
	if serviceName == "CANOES":
		# cmd = "docker run -dti --name="+containerName+" -v "+run+":"+run+" -v "+genome+":"+genome+" -v "+annotsvServer+":"+annotsvContainer+" "+serviceDockerImage+" /bin/bash /app/bin/canoes run -r "+run+" -s "+samplesheet+" -l "+bed+" -g "+genome+"/hg19.fa -o "+run
		cmd = "docker run --rm -dti --name="+containerName+" -v "+run+":"+run+" -v "+genome+":"+genome+" -v "+annotsvServer+":"+annotsvContainer+" "+serviceDockerImage+" /bin/bash /app/bin/canoes run -r "+run+" -s "+samplesheet+" -l "+bed+" -g "+genome+"/hg19.fa -o "+run
	subprocess.call(cmd, shell = True)
	createRunningFile(run, serviceName)
	createContainerFile(config, run, serviceName)
