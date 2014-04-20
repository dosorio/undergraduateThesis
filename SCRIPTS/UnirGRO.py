# Copyright 20013 © All rights reserved
# Paola Rondón-Villareal
# UnirGRO.py une los .gro de los péptidos y las membranas en un unico archivo.


# -*- coding: utf-8 -*-

fr=open('system.gro','r')

#Lectura para determinar el número de atomos

cont=-2
for line in fr:
    if cont==-2:
        line=next(fr)
        num1=int(line)
    
    cont=cont+1
    
    if cont==num1:
        line=next(fr)
        line=next(fr)
        num2=int(line)
        
fr.close()

#Lectura para crear el archivo que se necesita
fr=open('system.gro','r')
fw=open('system2.gro','w')

fw.write('Fixed')
fw.write("\n")

cont=-1
for line in fr:
    if cont==-1:
        line=next(fr)
        line=next(fr)
        numero=num1+num2
        print numero
        fw.write(str(numero))
        fw.write("\n")
        
        
    cont=cont+1
    
    if cont==num1:
        line=next(fr)
        line=next(fr)
        line=next(fr)
        cont=cont+1
    
    fw.write(line)
    #fw.write("\n")
        
fr.close()

