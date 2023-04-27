#!/usr/bin/env python3
#updated: 24 Apr 15
# Updated by Sian Bray on 12th November 2022

# Hash out lines 66 and 76 for absolute allele frequency difference
# Alternatively hash out lines 67 and 77 for raw allele frequency difference

import os, sys, subprocess

cohort1_name = sys.argv[1]
cohort1_n = sys.argv[2]

cohort2_name = sys.argv[3]
cohort2_n = sys.argv[4]

orient_genefile = sys.argv[5]

candidate_genefile = open(sys.argv[6],'r')

winsize = sys.argv[7]

count = 0

AF_list = []
for file in os.listdir():
    if file.endswith('.table'):
        AF_list.append(file)

# print('AF List:')
# print(AF_list)
# print('\n')

for line in candidate_genefile:
    data = line.split()
    locus = data[0]
    print('Processing locus '+locus)
    genename = data[1]
    scaffold = data[2]+'.table'
    genestart = int(data[3])
    geneend = int(data[4])
    genemiddle = int(genestart+((geneend-genestart)/2))

    count += 1

    for AF in AF_list:
        # print(scaffold)
        # print(AF)
        if scaffold in AF:
            infile = open(AF,'r')
            Rfile = open(locus+'_'+winsize+'.R','w')
            Rfile.write('s=read.table("'+AF+'",header=T)\n'+
                        'f=read.table("'+orient_genefile+'",header=F)\n\n'+
                        'locus="'+locus+'"\n'+
                        'genename="'+genename+'"\n'+
                        'scaffold="'+data[2]+'"\n'+
                        'genestart='+str(genestart)+'\n'+
                        'geneend='+str(geneend)+'\n'+
                        'genemiddle='+str(genemiddle)+'\n'+
                        'winsize='+winsize+'\n'+
                        'start=genemiddle-winsize\n'+
                        'end=genemiddle+winsize\n'+
                        's[is.na(s[,3]),3]=0\n'+ #if AC for Pop1 is NA change it to 0
                        's[is.na(s[,7]),7]=0\n'+ #if AC for Pop2 is NA change it to 0
                        's[,3]=s[,3]/'+cohort1_n+'\n'+ #Change AC to AF for Pop1
                        's[,7]=s[,7]/'+cohort2_n+'\n'+ #Change AC to AF for Pop2
                        #'s[,7]=s[,3]-s[,7]\n'+ # Allele frequency difference
                        's[,7]=abs(s[,3]-s[,7])\n'+ # Absolute allele frequency difference
                        'sr=s[s[,2]>=start & s[,2]<=end,]\n\n'+ #Select only sites that are larger than the start and smaller than the end
                        'chr=subset(f[,2],as.vector(f[,1])==locus)\n'+
                        'locstart=subset(f[,3],as.vector(f[,1])==locus)\n'+
                        'locend=subset(f[,4],as.vector(f[,1])==locus)\n'+
                        'genes=subset(f,(f[,2]==chr & f[,3]>=start & f[,3]<=end) | (f[,2]==chr & f[,4]>=start & f[,4]<=end))\n'+
                        'genes[,6]=1\n'+
                        'genes[as.vector(genes[,5])=="+",6]=2\n\n'+
                        'pdf(file=paste(locus,"_",winsize,".pdf",sep=""))\n'+
                        #'plot(sr[,2],sr[,7],bty="l",type="p",cex=1.0, col="black",pch=20,ylim=c(-1,1.2),ylab="Allele Frequency Difference ('+cohort1_name+' - '+cohort2_name+')",xlab="'+scaffold+' position (bp)")\n'+ # Allele frequency difference
                        'plot(sr[,2],sr[,7],bty="l",type="p",cex=1.0, col="black",pch=20,ylim=c(0,1.2),ylab="Absolute allele Frequency Difference.",xlab="'+scaffold+' position (bp)")\n'+ # Absolute allele frequency difference
                        'for(i in 1:nrow(genes)) {\n'+
                        '\tarrows(x0=genes[i,3],y0=1.1,x1=genes[i,4],code=genes[i,6],length=0.1,y1=1.1,col="grey",lwd=3)\n'+
                        '\t}\n'+
                        'text(x=genemiddle,y=1.20,labels=locus,col="red")\n'+
                        'dev.off()\n\n'+
                        'sum=list(locus,genename,scaffold,genestart,geneend,start,end,length(sr[,2]),max(sr[,7]),min(sr[,7]))\n'+
                        'write.table(sum,file=paste(locus,".sum",sep=""),sep="\t",row.names=F,col.names=F,quote=F)\n')
            Rfile.close()

            cmd = ('Rscript '+locus+'_'+winsize+'.R')
            p = subprocess.Popen(cmd, shell=True)
            sts = os.waitpid(p.pid, 0)[1]

            infile.close()

print('\nProcessed '+str(count)+' genes')

header = open('header','w')
header.write('Locus\tGeneName\tScaffold\tGene_start\tGene_end\tWindow_start\tWindow_end\tNum_SNPs\tMax_diff\tMin_diff\n')
header.close()

cmd = ('cat header *.sum > '+sys.argv[6].replace('.txt','.out'))
p = subprocess.Popen(cmd, shell=True)
sts = os.waitpid(p.pid, 0)[1]

cmd = ('rm header *.sum')
p = subprocess.Popen(cmd, shell=True)
sts = os.waitpid(p.pid, 0)[1]

print('\n\nFinished!!\n')
