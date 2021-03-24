#!/usr/bin/python

import sys,os
from optparse import OptionParser
from optparse import OptionGroup

import numpy as np

titv=15.0;

parser = OptionParser(usage="usage: %prog [options] [mutation.config] [reference.fa]");
parser.add_option("", "--titv",  dest="titv", help="Transition/transversion ratio default:("+str(titv)+")",default=titv,   type="float");
#group = OptionGroup(parser, "Group1","First Group of parameters");
#group.add_option("-o", "--outfile", dest="outfile", help="output");

if(len(sys.argv) == 1):
  parser.print_help()
  sys.exit(1)

  
#parser.add_option_group(group);
(options, args) = parser.parse_args();

if( len(args) < 2 ):
    sys.stderr.write("\nneed a least 2 arguments, please use -h to see options.\n");
    sys.exit(1);


conffile      = str(sys.argv[1]);
reffile       = str(sys.argv[2]);



filefa = open(reffile, "r");
linec=0;
header="";
genome="";
#read genome
for linefa in filefa:
  linefa = linefa.strip();

  if(len(linefa)==0):
    continue;
  #print(linefa);
  if(linec==0):
    header=linefa;
  else:
    genome=genome+linefa;

    
  linec+=1;
filefa.close();


lengenome=len(genome);


filecf = open(conffile, "r");
linec=0;
#read genome
genomicRange = [];


for linecf in filecf:

  print(linecf);
  #removing #
  linecf = linecf.split("#", 1)[0]
  print(linecf);
  linecf = linecf.strip();

  if(len(linecf)==0):
    continue;
  fields = linecf.split();
  
  if(len(fields)!=11):
    sys.stderr.write("\nThe line ->"+str(linecf)+"<- does not have 10 columns");#, please check the expected format using -h\n");
    sys.exit(1);

  rangechr = fields[0].split("-");
  if(fields[1] != "sn"):
    sys.stderr.write("field#2 should be sn");#, please check the expected format using -h\n");
    sys.exit(1);
  if(fields[3] != "in"):
    sys.stderr.write("field#4 should be in");#, please check the expected format using -h\n");
    sys.exit(1);
  if(fields[5] != "inl"):
    sys.stderr.write("field#6 should be inl");#, please check the expected format using -h\n");
    sys.exit(1);
  if(fields[7] != "de"):
    sys.stderr.write("field#8 should be de");#, please check the expected format using -h\n");
    sys.exit(1);
  if(fields[9] != "del"):
    sys.stderr.write("field#10 should be del");#, please check the expected format using -h\n");
    sys.exit(1);
    


  genomicRange.append( {
    "st":  int(rangechr[0]),
    "end":    int(rangechr[1]),
    "sn":     float(fields[2]),
    "in":     float(fields[4]),
    "inl":    float(fields[6]),
    "de":     float(fields[8]),
    "del":    float(fields[10])
  } );

  
  linec+=1;
filefa.close();


if(genomicRange[0]["st"] != 1):
    sys.stderr.write("first coordinate should be 1");#, please check the expected format using -h\n");
    sys.exit(1);

if(genomicRange[len(genomicRange)-1]["end"] != len(genome)):
    sys.stderr.write("last coordinate should be "+str(len(genome)));#, please check the expected format using -h\n");
    sys.exit(1);
  
for i in range(0,len(genomicRange)-1):
  if( genomicRange[i]["end"] != (genomicRange[i+1]["st"]-1) ):
    sys.stderr.write("coordinate "+str(i)+" should be the end of the other "+str( genomicRange[i]["st"] ) +" "+str(genomicRange[i+1]["end"]-1) );
    sys.exit(1);
    

