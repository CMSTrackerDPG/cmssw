#!/usr/bin/env python
import Config
import ConfigParser
import ROOT
import calendar
import commands
import datetime
import glob
import json
import numpy
import optparse
import os
import re
import string
import subprocess
import sys
import time
import urllib

from datetime import datetime
import CondCore.Utilities.conddblib as conddb

def log(s):
    print time.strftime("%Y-%m-%d-%H-%M-%S: ") + s

def getRunEndTime(con,run):

    con = conddb.connect(url = conddb.make_url())
    session = con.session()
    RunInfo = session.get_dbtype(conddb.RunInfo)

    bestRun = session.query(RunInfo.run_number,RunInfo.start_time, RunInfo.end_time).filter(RunInfo.run_number == run).first()
    if bestRun is None:
        raise Exception("Run %s can't be matched with an existing run in the database." %run)

    stop = bestRun[2]
    bestRunStopTime  = calendar.timegm( bestRun[2].utctimetuple() ) << 32

    #print "run stop time: ",stop,"(",bestRunStopTime,")"

    return bestRunStopTime

def getRunStartTime(con,run):

    con = conddb.connect(url = conddb.make_url())
    session = con.session()
    RunInfo = session.get_dbtype(conddb.RunInfo)

    bestRun = session.query(RunInfo.run_number,RunInfo.start_time, RunInfo.end_time).filter(RunInfo.run_number == run).first()
    if bestRun is None:
        raise Exception("Run %s can't be matched with an existing run in the database." %run)

    start= bestRun[1]
    bestRunStartTime = calendar.timegm( bestRun[1].utctimetuple() ) << 32

    #print "run start time:",start,"(",bestRunStartTime,")"
    return bestRunStartTime

def unpackLumiId(since):
    kLowMask = 0XFFFFFFFF
    run  = since >> 32
    lumi = since & kLowMask
    return run, lumi

def _high(n):
    return int(n) >> 32

def _low(n):
    return int(n) & 0xffffffff

def since_to_date(since,isDateFormat=False):
    """ converts tag since value to miliseconds """
    dt = datetime.utcfromtimestamp(_high(since))
    epoch = datetime.utcfromtimestamp(0) # January 1, 1970
    ms = (dt - epoch).total_seconds()
    if(not isDateFormat):
        return str(datetime.utcfromtimestamp(ms))
    else:
        return datetime.utcfromtimestamp(ms)

def isRunOneDayOld(con,run):
    endOfRun = getRunEndTime(con,run)
    #now = datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")

    delta = datetime.utcnow()-since_to_date(endOfRun,True)
    log("time difference {}".format(delta))
    if( delta.days > 1):
        return True
    else:
        return False

def isRunCoveredForTag(con,run,tag_name,myRecord):

    session = con.session()
    IOV     = session.get_dbtype(conddb.IOV)
    TAG     = session.get_dbtype(conddb.Tag)
    GT      = session.get_dbtype(conddb.GlobalTag)
    GTMAP   = session.get_dbtype(conddb.GlobalTagMap)
    RUNINFO = session.get_dbtype(conddb.RunInfo)

    RefTagIOVs = session.query(IOV.since,IOV.payload_hash,IOV.insertion_time).filter(IOV.tag_name ==tag_name).all()
    TagInfo = session.query(TAG.synchronization,TAG.time_type).filter(TAG.name == tag_name).all()[0]

    if(TagInfo[1]=='Run'):
        log("Time type is: {}".format(TagInfo[1]))
        for i in RefTagIOVs:
            if ( str(i[0])==str(run) ):
                log("Found an IOV for run ={}".format(run))
                return True
    elif (TagInfo[1]=='Lumi'):
        log("Time type is {}".format(TagInfo[1]))
        for i in RefTagIOVs:
            if(unpackLumiId(i[0])[0]>run):
                log("Found an IOV ({}) > than the run={}".format(unpackLumiId(i[0]),run))
                return True
    elif (TagInfo[1]=='Time'):
        end_time   = getRunEndTime(con,run)
        log("Time type is: {}".format(TagInfo[1]))
        for i in RefTagIOVs:
            newdate1 = time.strptime(since_to_date(i[0]),"%Y-%m-%d %H:%M:%S")
            newdate2 = time.strptime(since_to_date(end_time),"%Y-%m-%d %H:%M:%S")
            if(newdate1 > newdate2):
                log("Found an IOV ({}) > that the run endtime={}".format(since_to_date(i[0]),since_to_date(end_time)))
                return True
    else:
        log("Should never happen!")


def isStreamDone(con,run):
  if(isRunOneDayOld(con,run)):
    print "WARNING, FORCE PROCESSING %s AS IT IS MORE THAN 1D OLD"%run
    return True

  out = subprocess.check_output(["curl", "-k", "-s", "https://cmsweb.cern.ch/t0wmadatasvc/prod/run_stream_done?run={}&stream=Express".format(run)])
  m = re.match('{"result": \[\n (.*)\]}\n', out)
  if m:
    return m.group(1) == "true"
  else:
    print "Could not get correct info for run", run
    return False

def dasQuery(query,config):
  cmd = config.dasClient+" --limit=9999 --query=\"%s\""%query
  print cmd
  output = commands.getstatusoutput(config.initEnv+cmd)
  if output[0]!=0:
    print "DAS CLIENT ERROR %s"%output[0]
    print output[1]
    sys.exit()
  return(output[1].splitlines())

def checkDatasetStructure(dataset,silent=False):
   goodDataset=True
   d = dataset.split("/")
   if not len(d) == 4:
      if not silent:
         print "Bad dataset. Expecting 4 '/'"
      goodDataset=False
      return False
   if not d[0]=='':
      if not silent:
         print "Bad dataset. Expecting nothing before first '/'"
      goodDataset=False
   if not len(d[1])>0 or not len(d[2]) > 0 or not len(d[3]) > 0:
      if not silent:
         print "Bad dataset. Expecting text between '/'"
      goodDataset=False
   if os.path.isdir(dataset):
      if not silent:
         print "Bad dataset. Can't be an existing directory"
      goodDataset=False
   return goodDataset

def getDatasetFromPattern(pattern,conf):
   if not checkDatasetStructure(pattern):
      print "FATAL ERROR, bad dataset pattern"
      return([])
   result = dasQuery("dataset dataset=%s"%pattern,conf)
   datasets = []
   for line in result:
      print line
      if checkDatasetStructure(line,silent=False):
         datasets.append(line)
   return datasets

def getRunsFromDataset(dataset,conf):
   if not checkDatasetStructure(dataset):
      print "FATAL ERROR, bad dataset pattern"
      return([])
   result = dasQuery("run dataset=%s"%dataset,conf)
   runs=[]
   for line in result:
      if line.isdigit:
         if len(line)==6: #We want the run number to be btw 100.000 and 999.999
            runs.append(int(line))
   runs.sort()
   return runs

def getNumberOfEvents(run,dataset,conf):
   if not int(run) > 99999 or not int(run)<1000000:
      print "Invalid run number"
      return 0
   if not checkDatasetStructure(dataset):
      print "Invalid dataset"
      return 0



   NEventsDasOut = dasQuery("summary run=%s dataset=%s |grep summary.nevents"%(run,dataset),conf)[-1].replace(" ","")
   if not NEventsDasOut.isdigit():
      print "Invalid number of events:"
      print "__%s__"%NEventsDasOut
      return 0
   else:
      return int(NEventsDasOut)

def getNumberOfFiles (run,dataset,conf):
   if not int(run) > 99999 or not int(run)<1000000:
      print "Invalid run number"
      return 0
   if not checkDatasetStructure(dataset):
      print "Invalid dataset"
      return 0
   NFilesDasOut = dasQuery('summary dataset=%s run=%s | grep summary.nfiles'%(dataset,run),conf)[-1].replace(" ","")
   if not NFilesDasOut.isdigit():
      print "Invalid number of files."
      return 0
   else :
      return int(NFilesDasOut)


def reSubmitJob(run, dataset, conf, first, last):
   print "Re-submitting jobs for run = %s, dataset = %s"%(run, dataset)



def submitJobs(run, dataset, nFiles, conf):
   print "Submitting jobs for run = %s, dataset = %s"%(run, dataset)

   #GET THE LIST OF FILE FROM THE DATABASE
   files = ''
   if not checkDatasetStructure(dataset,conf):
      print "FATAL ERROR, bad dataset"
      return -1
   if not run > 99999 or not run<1000000:
      print "FATAL ERROR, bad run number"
      return -1
   filesList = dasQuery('file dataset=%s run=%s'%(dataset,run),conf)
   filesInJob = 0
   firstFile = 0
   for f in filesList:
      if(not f.startswith('/store')):continue
      if filesInJob<conf.nFilesPerJob:
         files+="'"+f+"',"
         filesInJob+=1
      else:
         sendJob(dataset,run,files,conf,firstFile)
         firstFile = firstFile+filesInJob
         files="'"+f+"',"
         filesInJob=1
   sendJob(dataset,run,files,conf,firstFile)

def sendJob(dataset,run,files,conf,firstFile):
   cmd = "python %s/submitCalibTree/runJob.py -f %s --firstFile %s -d %s -r %s "%(conf.RUNDIR, files,firstFile,dataset,run)
   if conf.AAG:
      cmd+=" -a "
   bsub = 'bsub -q 2nd -J calibTree_' + str(run) + '_' + str(firstFile)+ '_' + '_%s'%("Aag" if conf.AAG else 'Std')+' -R "type == SLC6_64 && pool > 30000" ' + ' "'+cmd+'"'
   conf.launchedRuns.append([run,firstFile])
   if conf.submit:
      os.system(bsub)
   else:
      print cmd + " --stageout False"

def generateJobs(conf):
   print "Gathering jobs to launch."
   print conf
   lastRunProcessed = conf.firstRun
   datasets = getDatasetFromPattern(conf.datasetPat,conf)

   import CondCore.Utilities.conddblib as conddb
   con = conddb.connect(url = conddb.make_url("pro"))

   for d in datasets:
      datasetRuns = getRunsFromDataset(d,conf)
      print datasetRuns
      for r in datasetRuns:
         if int(r) > conf.firstRun and int(r)<conf.lastRun:
            print "Checking run %s"%r

            hasDetVOff = isRunCoveredForTag(con,r,"SiStripDetVOff_13hourDelay_v1_Validation","SiStripDetVOffRcd")
            if(hasDetVOff):
              log("run {} has DetVOff coverage".format(r))
            else:
              log("run {} is NOT covered by DetVOff".format(r))
              break

            if not isStreamDone(con,r):
              print "Stream not processed yet... break."
              break

            n=getNumberOfEvents(r,d,conf)
            if n < 250:
               print "Skipped. (%s evt)"%n
            else:
               nFiles = getNumberOfFiles(r,d,conf)
               if nFiles > 0:
                  print "Will be processed ! (%s evt, %s files)"%(n,nFiles)
                  if r > lastRunProcessed:
                     lastRunProcessed = r
                  submitJobs(r,d,nFiles,conf)
               else:
                  print "Skipped. (%s evt,%s files)"%(n,nFiles)
         else:
            for failled in conf.relaunchList:
               if int(failled[0]) == int(r):
                  print "Relaunching job %s "% failled
                  if len(failled)==3:
                     reSubmitJob(int(failled[0]),d,conf,failled[1],failled[2])
                  else:
                     submitJobs(int(failled[0]),d,25,conf)
   return lastRunProcessed


def cleanUp():
   os.system('rm core.*')


def checkCorrupted(lastGood, config):
   calibTreeList = ""
   print("Get the list of calibTree from" + config.CASTORDIR + ")")
   calibTreeInfo = commands.getstatusoutput(config.eosLs +config.CASTORDIR)[1].split('\n');
   NTotalEvents = 0;
   runList = []

   for info in calibTreeInfo:
      subParts = info.split();
      if(len(subParts)<4):  continue

      runList.append(subParts[-1].replace("calibTree_","").replace(".root","").split("_"))
   print runList
   datasets = getDatasetFromPattern(config.datasetPat,config)
   for d in datasets:
      datasetRuns = getRunsFromDataset(d,config)
      print datasetRuns
      for r in datasetRuns:
         if int(r) > lastGood:
            print "Checking run %s"%r
            n=getNumberOfEvents(r,d,config)
            if n < 250:
               print "Skipped. (%s evt)"%n
            else:
               nFiles = getNumberOfFiles(r,d,config)
               if nFiles < 25:
                  print "Found run %s ? %s"%(r,[str(r)] in runList)
               else:
                  x=25
                  while x<nFiles:
                     print "Found run %s , %s ? %s "%(r,x, [str(r),str(x)] in runList)
                     x+=25





   #for line in runList:
   #   print line
      #file = "root://eoscms/"+config.CASTORDIR+"/"+subParts[-1]
      #print("Checking " + file)
      #results = commands.getstatusoutput(config.initEnv+'root -l -b -q ' + file)
      #if(len(results[1].splitlines())>3):
      #   print(results[1]);
      #   print("add " + str(run) + " to the list of failled runs")
#         os.system('echo ' + str(run) + ' >> FailledRun%s.txt'%('_Aag' if AAG else ''))



if __name__ == "__main__":
   print "DEBUG DEBUG DEBUG"
   c = Config.configuration(False)
   c.runNumber = 1
   c.firstRun  = 274500
   c.lastRun   = 275000
   c.debug     = True
   #generateJobs(c)
   checkCorrupted(0,c)

if False:
#elif(checkCorrupted):
   #### FIND ALL CORRUPTED FILES ON CASTOR AND MARK THEM AS FAILLED RUN

   calibTreeList = ""
   print("Get the list of calibTree from" + CASTORDIR + ")")
   calibTreeInfo = commands.getstatusoutput(initEnv+"eos ls " + CASTORDIR)[1].split('\n');
   NTotalEvents = 0;
   run = 0
   for info in calibTreeInfo:
      subParts = info.split();
      if(len(subParts)<4):continue

      run = int(subParts[4].replace("/calibTree_","").replace(".root","").replace(CASTORDIR,""))
      file = "root://eoscms//eos/cms"+subParts[4]
      print("Checking " + file)
      results = commands.getstatusoutput(initEnv+'root -l -b -q ' + file)
      if(len(results[1].splitlines())>3):
         print(results[1]);
         print("add " + str(run) + " to the list of failled runs")
         os.system('echo ' + str(run) + ' >> FailledRun%s.txt'%('_Aag' if AAG else ''))

   #### If mode = All, relaunch with mode = Aag
   if opt.datasetType.lower()=="all":
      system("cd "+RUNDIR+"; python SubmitJobs.py -c -d Aag")

#else:
   #### UNKNOWN CASE
#   print "unknown argument: make sure you know what you are doing?"
