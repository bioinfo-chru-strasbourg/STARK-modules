#!/usr/bin/env python

from flask import Flask, request, jsonify
import sys
import optparse
import time
import os
import subprocess
import urllib
import json
import re
import string
import random



# VARIABLES
#############


# Docker STARK image

if 'DOCKER_STARK_IMAGE' in os.environ:
    docker_stark=os.environ['DOCKER_STARK_IMAGE']
else:
    docker_stark="stark"

# Task Spooler variables

if 'TS' in os.environ:
    ts=os.environ['TS']
else:
    ts=""

if 'TS_SOCKET' in os.environ:
    ts_socket=os.environ['TS_SOCKET']
else:
    ts_socket="/ts-tmp/STARK"

if 'TS_SAVELIST' in os.environ:
    ts_savelist=os.environ['TS_SAVELIST']
else:
    ts_savelist="/ts-tmp"

if 'SHELL' in os.environ:
    shell=os.environ['SHELL']
else:
    shell="/bin/ash"


#ts_env=" TS_SOCKET=" + ts_socket + " TS_SAVELIST=" + ts_savelist + " "
ts_env=" TS_SAVELIST=" + ts_savelist + " "


# DOCKER STARK CONTAINER MOUNT

if 'DOCKER_STARK_SERVICE_STARK_API_CONTAINER_MOUNT' in os.environ:
    docker_stark_container_mount=os.environ['DOCKER_STARK_SERVICE_STARK_API_CONTAINER_MOUNT']
else:
    docker_stark_container_mount=""


# STARK API LOG Folder 

if 'DOCKER_STARK_SERVICE_STARK_API_LOG_FOLDER' in os.environ:
    docker_stark_api_log_folder=os.environ['DOCKER_STARK_SERVICE_STARK_API_LOG_FOLDER']
else:
    docker_stark_api_log_folder="/STARK/services/stark/stark/api"

# STARK API RUNS Folder 

if 'DOCKER_STARK_SERVICE_STARK_API_RUNS_FOLDER' in os.environ:
    docker_stark_api_runs_folder=os.environ['DOCKER_STARK_SERVICE_STARK_API_RUNS_FOLDER']
else:
    docker_stark_api_runs_folder="/STARK/input/runs"


# HOME
home_header="************* \n"
home_header+="**<B>STARK API</B>** \n"
home_header+="************* \n"
home_header+="<A href='help'>Help</A> \n"
home_header+="<A href='queue?help&format=html'>List of task</A> \n"

# HELP
help_header="******** \n"
help_header+="**<B><A href='/help'>HELP</A></B>** \n"
help_header+="******** \n"

help_queue="**<B>queue</B>** \n"
help_queue+="/queue?list: list all task \n"
help_queue+="/queue?info=ID: task info file (default running/last task) \n"
help_queue+="/queue?log=ID: task log file (default running/last task) \n"
help_queue+="/queue?kill=ID: kill task (default running task) \n"
help_queue+="/queue?prioritize=ID: task prioritization (default last task) \n"
help_queue+="/queue?remove=ID: remove task (default last task) \n"
help_queue+="/swap=ID1-ID2: swap 2 tasks ID1 and ID2 \n"
help_queue+="Add '&format=html' to view in html\n"

help_analysis="**<B>analysis</B>** \n"
help_analysis+="/analysis: Launch STARK analysis with JSON parameters in POST mode\n"
help_analysis+="Command example with curl in POST method: curl -s -X POST -H 'Content-Type: application/json' -d '{\"run\":\"MY_RUN\"}' http://IP:PORT/analysis \n"
help_analysis+="Command example with URL in GET method: http://IP:PORT/analysis?json={\"run\":\"MY_RUN\"}  $API_URL \n"
help_analysis+="Run launch (run name detected if exists) JSON example: '{\"run\":\"MY_RUN\"}' \n"
help_analysis+="Help JSON example: '{\"help\":null}' \n"
help_analysis+="Application JSON example: '{\"applications_infos\":null}' \n"


# DEF

def randomStringDigits(stringLength=6):
    """Generate a random string of letters and digits """
    lettersAndDigits = string.ascii_letters + string.digits
    return ''.join(random.choice(lettersAndDigits) for i in range(stringLength))



app = Flask(__name__)

#app.url_map.strict_slashes = True

start = int(round(time.time()))


@app.route("/", methods = ['GET','POST']) #methods = ['GET', 'POST', 'DELETE'])
def default():
    ret='<pre>'+home_header+'\n'+help_header+'\n'+help_queue+'\n'+help_analysis+'</pre>'
    return ret, 200


@app.route("/help", methods = ['GET','POST']) #methods = ['GET', 'POST', 'DELETE'])
def help():
    ret='<pre>'+home_header+'\n'+help_header+'\n'+help_queue+'\n'+help_analysis+'</pre>'
    return ret, 200




@app.route("/analysis", methods = ['GET','POST']) #methods = ['GET', 'POST', 'DELETE'])
#def hello_world(json=None):
def stark_launch():

    # Example:
    # curl -s -X POST -H 'Content-Type: application/json' -d '{"run":"run_name"}' $LAUNCHER
    
    # Main param
    ##############

    # param
    input_method=""
    header=""

    # Input JSON data
    json_input = request.get_json(force=False)
    json_dump=""+str(json.dumps(json_input))

    # Test input data
    # data in post
    if json_dump != "null":
        json_input = request.get_json(force=True)
        json_dump=json.dumps(json_input)
        input_method="POST"
    # data in get
    elif 'json' in request.args:
        json_dump=request.args['json']
        json_input = json.loads(json_dump)
        input_method="GET"
        header=home_header+'\n'+'STARK ANALYSIS ID: '
        if 'format' in request.args:
            if request.args['format']=="html":
                header="<pre>\n"+header+"</pre>\n"
    # HELP
    elif 'help' in request.args:
        ret=home_header+'\n'+help_header+'\n'+help_analysis
        if 'format' in request.args:
            if request.args['format']=="html":
                ret="<pre>\n"+ret+"</pre>\n"
        return ret, 200
    # no data, HELP
    else:
        ret=home_header+'\n'+help_header+'\n'+help_analysis
        if 'format' in request.args:
            if request.args['format']=="html":
                ret="<pre>\n"+ret+"</pre>\n"
        return ret, 200


    # analysis ID & NAME default
    analysesID=randomStringDigits(12)
    analysesRUNNAME=analysesID
    analysesNAME=analysesID


    # Folders
    ###########

    # analysisID
    ##############

    # Find analysisID from RUN info in JSON

    runID=""
    runMD5=""
    if 'run' in json_input.keys():
        #print json_input['run']
        if len(json_input['run'].split(",")) == 1:
            
            # Run infos
            name=json_input['run'].split(":")[0]
            runID=os.path.basename(name)
            runFolder=""
            if os.path.isdir(name):
                runFolder=name
            # elif os.path.isdir(DOCKER_STARK_INNER_FOLDER_INPUT_RUNS+"/"+name):
            #     runFolder=DOCKER_STARK_INNER_FOLDER_INPUT_RUNS+"/"+name
            elif os.path.isdir(docker_stark_api_runs_folder+"/"+name):
                runFolder=docker_stark_api_runs_folder+"/"+name

            # Command to find MD5/sha1sum
            if runFolder!="":
                myCmd ="find " + runFolder + " -maxdepth 1 -type f -print0 | xargs -0 sha1sum | cut -b-40 | sha1sum | awk '{print $1}'"
            else:
                myCmd ="echo " + name + " | sha1sum | awk '{print $1}'"

            runMD5 = subprocess.check_output(myCmd, shell=True).strip()

    if runID != "" and runMD5 !="":
        #analysesNAME=runID
        analysesRUNNAME=runID
        analysesNAME="ID-" + runMD5 + "-NAME-" + runID


    # Find analysisID from ANALYSIS_NAME info in JSON

    analysis_name=""
    if 'analysis_name' in json_input.keys():
        #print json_input['run']
        analysis_name=json_input['analysis_name']

    if analysis_name != "":
        analysesRUNNAME=analysis_name
        analysesNAME=analysis_name


    # DOCKER parameters
    #####################

    # DOCKER MOUNT parameters
    
    docker_mount=docker_stark_container_mount


    # analysisIDNAME

    analysisIDNAME='STARK.' + analysesID + '.' + analysesNAME
    PostanalysisIDNAME='STARK-POSTANALYSIS.' + analysesID + '.' + analysesNAME


    # DOCKER name

    docker_name=" --name " + analysisIDNAME + " "


    # DOCKER parameters

    #docker_parameters=' --rm ' + docker_mount + ' ' + docker_env + ' ' + docker_name + ' '
    docker_parameters=' --rm ' + docker_mount + ' ' + docker_name + ' '


    # analysis FOLDERS et files
    
    #analysisFOLDER=DOCKER_STARK_INNER_FOLDER_SERVICES+"/"+DOCKER_STARK_SERVICES_FOLDER_STARK_API
    analysisFOLDER=docker_stark_api_log_folder
    analysisFILE=analysisFOLDER+"/"+analysisIDNAME+".json"
    analysisLOG=analysisFOLDER+"/"+analysisIDNAME+".log"
    analysisERR=analysisFOLDER+"/"+analysisIDNAME+".err"
    analysisINFO=analysisFOLDER+"/"+analysisIDNAME+".info"
    analysisOUTPUT=analysisFOLDER+"/"+analysisIDNAME+".output"
    analysisLOGERR_PARAM=' 1>' + analysisLOG + ' 2>' + analysisERR
    f = open(analysisFILE, "w")
    f.write(json_dump)
    f.close()


    # TASK SPOOLER

    if ts != "":
        ts_cmd=ts_env + ts + ' -L ' + analysisIDNAME
        #ts_cmd= ts + ' -L ' + analysisIDNAME
    else:
        ts_cmd=""

    
    # Prepare and Launch TS Docker run command
    
    myCmd = ts_cmd + ' docker run ' + docker_parameters + ' ' + docker_stark + ' --analysis_name=' + analysesRUNNAME + ' --analysis=' + analysisFILE + ' '#+ analysisLOGERR_PARAM
    #print(myCmd);
    getCmd = subprocess.check_output(myCmd, shell=True)
    STARKCmdID=getCmd.strip()

    # Retrive INFO from command
    # Example :
    # FIRST_TASKID=`ts ash -c "sleep 10; echo hi"`
    # ts ash -c "ts -w $FIRST_TASKID && echo there"
    if ts != "":
        # INFO
        #myPostCmd=ts_env + " " + ts + ' -L ' + PostanalysisIDNAME + ' ' + shell + ' -c "' + ts_env + ' && ' + ts + ' -i ' + STARKCmdID + ' > ' + analysisINFO + '"'
        #getPostCmd = subprocess.check_output(myPostCmd, shell=True);
        # OUTPUT
        #myPostCmd=ts_env + " " + ts + ' -L ' + PostanalysisIDNAME + ' ' + shell + ' -c "' + ts_env + ' && ' + ts + ' -c ' + STARKCmdID + ' > ' + analysisOUTPUT + '"'
        #getPostCmd = subprocess.check_output(myPostCmd, shell=True);
        # INFO && OUTPUT
        myPostCmd=ts_env + " " + ts + ' -L ' + PostanalysisIDNAME + ' ' + shell + ' -c "' + ts_env + ' && ' + ts + ' -i ' + STARKCmdID + ' > ' + analysisINFO + ' && ' + ts + ' -c ' + STARKCmdID + ' > ' + analysisOUTPUT + '"'
        getPostCmd = subprocess.check_output(myPostCmd, shell=True)


    # Notification or something else
    # command to launch : ( ts -w ; xmessage Finished! ) &

    # Return

    if getCmd:
        return header+analysisIDNAME, 200
    else:
        return "KO", 400



@app.route("/queue", methods = ['GET','POST']) #methods = ['GET', 'POST', 'DELETE'])
def queue():

    # params
    header=""
    format=""

    # queue action
    if 'log' in request.args:
        id=request.args['log']
        myCmd = 'cat $(ts -o ' + str(id) + ')'
    elif 'kill' in request.args:
        id=request.args['kill']
        myCmd = 'ts -k ' + str(id) + ' && echo "#[INFO] task #' + id + ' killed" ' + ' &&  ts -l'
    elif 'info' in request.args:
        id=request.args['info']
        myCmd = 'ts -i ' + str(id)
    elif 'prioritize' in request.args:
        id=request.args['prioritize']
        myCmd = 'ts -u ' + str(id) + ' && echo "#[INFO] task #' + id + ' prioritized" ' + ' &&  ts -l'
    elif 'remove' in request.args:
        id=request.args['remove']
        myCmd = 'ts -r ' + str(id) + ' && echo "#[INFO] task #' + id + ' removed"'  + ' &&  ts -l'
    elif 'swap' in request.args:
        id=request.args['swap']
        myCmd = 'ts -U ' + str(id) + ' && echo "#[INFO] tasks #' + id + ' swapped" ' + ' &&  ts -l'
    elif 'list' in request.args:
        myCmd = 'ts -l'
    elif 'help' in request.args:
        myCmd = 'ts -l'
        format="html"
        header=help_header+'\n'+help_queue+'\n'
    else:
        myCmd = 'ts -l'
        format="html"


    #print myCmd

    # format
    if 'format' in request.args:
        format=request.args['format']
    #else:
    #    format=""

    # command
    try:
        getCmd = subprocess.check_output(myCmd, shell=True)
        ret_val=200
    except:
        getCmd = "#[ERROR] Command error"
        ret_val=500


    # output
    if format == "html":
        ret="<pre>"+home_header+'\n'+header+getCmd+"</pre>"
    else:
        ret=getCmd

    return ret, ret_val


if __name__ == '__main__':
    parser = optparse.OptionParser(usage="python simpleapp.py -p ")
    parser.add_option('-p', '--port', action='store', dest='port', help='The port to listen on.')
    (args, _) = parser.parse_args()
    if args.port == None:
        print("Missing required argument: -p/--port")
        sys.exit(1)
    app.run(host='0.0.0.0', port=int(args.port), debug=False)
