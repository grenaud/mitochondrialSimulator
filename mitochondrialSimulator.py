#!/usr/bin/python

import sys,os
from optparse import OptionParser
from optparse import OptionGroup

import numpy as np
import random
import re

titv=15.0;
basefreq=[0.0,0.0,0.0,0.0];
cumbasefreq=[0.0,0.0,0.0,0.0];
snpcount=0;
inscount=0
delcount=0

def ranbombase():
  p =   random.random();
  if( (0 <= p)              and (p < cumbasefreq[0]) ):
    return 'A';
  if( (cumbasefreq[0] <= p) and (p < cumbasefreq[1]) ):
    return 'C';
  if( (cumbasefreq[1] <= p) and (p < cumbasefreq[2]) ):
    return 'G';
  if( (cumbasefreq[2] <= p) and (p < cumbasefreq[3]) ):
    return 'T';
  sys.stderr.write("Could not generate base error.\n");
  sys.exit(1);

def ranbombases(n):
  strtoreturn="";
  for i in range(0,n):
    strtoreturn+=ranbombase();
  return strtoreturn;

def transition(b):
  if(b == 'A'):
    return 'G';
  if(b == 'G'):
    return 'A';
  if(b == 'C'):
    return 'T';
  if(b == 'T'):
    return 'C';
  if(b == 'N'):
    return 'N';

  sys.stderr.write("transition error.\n");
  sys.exit(1);

def transversion(b):
  if(b == 'A'):
    #return 'G';
    return ['C','T'][ random.randint(0,1)];
  if(b == 'G'):
    #return 'A';
    return ['C','T'][ random.randint(0,1)];
  if(b == 'C'):
    #return 'T';
    return ['A','G'][ random.randint(0,1)];
  if(b == 'T'):
    #return 'C';
    return ['A','G'][ random.randint(0,1)];
  if(b == 'N'):
    return 'N';

  sys.stderr.write("transversion error.\n");
  sys.exit(1);

  

def mutate(seq,sn,ins,insl,des,desl):
  global snpcount;
  global inscount;
  global delcount;
  newseq="";
  #print(seq);
  myiter = iter(range(0,len(seq)));

  for i in myiter:
    pmut =   random.random();
    #print(pmut);
    if(pmut < sn):
      #print(pmut);
      #print(sn);
      ptitv =   random.random();
      snpcount+=1;
      if(ptitv < (float(titv)  / (float(titv) + float(1.0) ))):#is a transition
         newseq+=transition(seq[i]);
      else:
         newseq+=transversion(seq[i]);               
    else:    
      if( ( (sn)     <= pmut ) and (pmut < (sn+ins)     ) ):
        newseq+=seq[i];

        #print("ins");
        lins = np.random.zipf(insl);
        sins = ranbombases( lins );
        #sys.stderr.write("ins:"+sins+"\n");

        newseq  += sins;
        inscount+=1;
      else:
        if( ( (sn+ins) <= pmut ) and (pmut < (sn+ins+des) ) ):
          #print("del");
          delcount+=1;
          ldel = np.random.zipf(desl);
          if(ldel>=len(seq)):
            sys.stderr.write("A deletion has deleted an entire segment.\n");
            sys.exit(1);

                
          if(i == ( len(seq)-ldel)):#we removed everything to the last char
            return newseq;
                
          if(i < ( len(seq)-ldel)):#we have characters to delete
            for j in range(0,ldel-1):
                next(myiter, None)
            continue;

          #we need to delete chars from newseq
          lcharleft = len(seq)-i;
          lchartodel = ldel - lcharleft;
          newseq=newseq[:-lchartodel]
          #we do not have any more chars, we break
          break;

        else:
          newseq+=seq[i];#nothing happened
  return newseq;


         
  
      
outprefix="out";
outname="sim";
gen=1;
parser = OptionParser(usage="usage: %prog [options] [mutation.config] [reference.fa]");
parser.add_option("", "--titv",  dest="titv", help="Transition/transversion ratio default:("+str(titv)+")",default=titv,   type="float");
parser.add_option("-o", "--out",  dest="outprefix", help="Out prefix (default: "+str(outprefix)+")",default=outprefix,   type="string");
parser.add_option("", "--name",   dest="outname", help="Name of sequence (default:"+str(outname)+")",default=outname,   type="string");
parser.add_option("-g", "--gen",   dest="gen", help="Number of generations to simulate (default:"+str(gen)+")",default=gen,   type="int");

#group = OptionGroup(parser, "Group1","First Group of parameters");
#group.add_option("-o", "--outfile", dest="outfile", help="output");

#for i  in range(0,1000):
#  print(np.random.zipf(5));
  
if(len(sys.argv) == 1):
  parser.print_help()
  sys.exit(1)

  
#parser.add_option_group(group);
(options, args) = parser.parse_args();

if( len(args) < 2 ):
    sys.stderr.write("\nneed a least 2 arguments, please use -h to see options.\n");
    sys.exit(1);


conffile      = str(args[0]);
reffile       = str(args[1]);



filefa = open(reffile, "r");
linec=0;
header="";
genome="";
linelen=0;
#read genome
for linefa in filefa:
  linefa = linefa.strip();

  if(len(linefa)==0):
    continue;
  #print(linefa);
  if(linec==0):
    header=linefa;
  else:
    genome=genome+linefa.upper();
  if(linec == 1):
    linelen=len(linefa);
    
  linec+=1;
filefa.close();


lengenome=len(genome);
for i in range(lengenome):
  if(genome[i] == 'A'):
    basefreq[0]+=1;
  if(genome[i] == 'C'):
    basefreq[1]+=1;
  if(genome[i] == 'G'):
    basefreq[2]+=1;
  if(genome[i] == 'T'):
    basefreq[3]+=1;


for i in range(0,4):
  basefreq[i]=basefreq[i]/float(lengenome);
  #print(basefreq[i]);

cumbasefreq[0]=basefreq[0]

for i in range(1,4):
  cumbasefreq[i]=cumbasefreq[i-1]+basefreq[i];
cumbasefreq[3]=1;


filecf = open(conffile, "r");
linec=0;
#read genome
genomicRange = [];


for linecf in filecf:

  #print(linecf);
  #removing #
  linecf = linecf.split("#", 1)[0]
  #print(linecf);
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
    

genomicSegments = [];

for i in range(0,len(genomicRange)):
  stc=genomicRange[i]["st"]-1;
  enc=genomicRange[i]["end"];            
  genomicSegments.append( genome[stc:enc]);
  #print( str(stc)+"\t"+str(enc));
  #print( genome[stc:enc]);

sanitytestg="";
for i in range(0,len(genomicSegments)):
  sanitytestg=sanitytestg+genomicSegments[i];


if(sanitytestg != genome):
    sys.stderr.write("The segmented genomes does not equal to the original one, probably due to some gap in the genomic ranges" );
    sys.exit(1);


for g in range(1,options.gen+1):
  genometowrite="";
  sys.stderr.write("Generation#"+str(g)+"\n" );
  snpcountT=0;
  inscountT=0
  delcountT=0
  
  for i in range(0,len(genomicRange)):
    snpcount=0;
    inscount=0
    delcount=0

    ng= mutate( genomicSegments[i], genomicRange[i]["sn"], genomicRange[i]["in"], genomicRange[i]["inl"], genomicRange[i]["de"], genomicRange[i]["del"]);
    genometowrite = genometowrite + ng;
    genomicSegments[i]=ng;
    snpcountT+=snpcount;
    inscountT+=inscount;
    delcountT+=delcount;
    sys.stderr.write(str(genomicRange[i]["st"])+"-"+str(genomicRange[i]["end"])+" "+str(genomicRange[i]["end"]-genomicRange[i]["st"])+"bp SNPS:"+str(snpcount)+" ins:"+str(inscount)+" del:"+str(delcount)+"\n" );
  sys.stderr.write("TOTAL SNPS:"+str(snpcountT)+" ins:"+str(inscountT)+" del:"+str(delcountT)+"\n" );
  genometowrite=re.sub("(.{"+str(linelen)+"})", "\\1\n", genometowrite, 0, re.DOTALL);
  
  outfastafp = open(options.outprefix+"_"+str(g)+".fa", "w");
  outfastafp.write(">"+options.outname+" "+str(g)+"\n");
  outfastafp.write(genometowrite+"\n");
  outfastafp.close();

sys.stderr.write("Done\n" );
