# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 10:45:46 2014

"""


import glob
import os
os.chdir("./")
for file in glob.glob("L-*.pdb"):
    inFile=open(file, 'r')
    outFile=open(file+'-A','w')
    
    for line in iter(inFile):
        
        if line.find('ATOM')>=0:
            name=line[12:16]
            temp=name[0]
            
            if temp.isdigit():
                newLine=line[0:12]+name[1:len(name)]+'\t' +' ' +line[17:len(line)]
                
                outFile.write(newLine)
                
            else:
                outFile.write(line)
        else:
            outFile.write(line)           
        
    inFile.close()
    outFile.close()
            
