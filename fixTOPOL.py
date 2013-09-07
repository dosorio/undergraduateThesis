# -*- coding: utf-8 -*-

fr=open('topol.top','r')
fw=open('topol2.top','w')

#Lectura del archivo linea por linea
for line in fr:
    newline=line
    
    if line.find("gromos53a6.ff/")>=0:
        #The line needs to be fixed
        #First part of the line
        ini=0
        fin=line.find("gromos53a6.ff/")
        sub1=line[ini:fin]
        #Second part of the line
        ini=fin+14
        fin=len(line)
        sub2=line[ini:fin]
        #newline
        newline=sub1+ "gromos53a6_lipid.ff/" + sub2
        print line
    
    fw.write(newline)
    
    if line.find("Include water topology")>=0:
        line=next(fr)
        #First part of the line
        ini=0
        fin=line.find("gromos53a6.ff/")
        sub1=line[ini:fin]
        #Second part of the line
        ini=fin+14
        fin=len(line)
        sub2=line[ini:fin]
        #newline
        newline=sub1+ "gromos53a6_lipid.ff/" + sub2
        fw.write(newline)
        fw.write("; Include POPC topology")
        fw.write("\n")
        fw.write("#include \"POPC.itp\"")
        fw.write("\n")
        
    if line.find("[ molecules ]")>=0:
        line=next(fr)
        fw.write(line)
        line=next(fr)
        fw.write(line)
        line="POPC	128"
        fw.write(line)
        
    
        
fr.close()
fw.close()

