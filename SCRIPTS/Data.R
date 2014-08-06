#Paquetes necesarios
require(Peptides)
require(Rpdb)

#Lugar de los archivos
setwd("    ")

#CALCULOS
#Distancia
setwd("/Users/Daniel/Desktop/BioInformatics_Data/")
archivos<-list.files("DIST/DIST_POPG/")
setwd("DIST/DIST_POPG/")
for (i in (1:60)){
  pep<-grep(as.vector(POPG$ID)[i],archivos)[1]
  POPG$minD[i]<-round(min(abs(read.xvg(archivos[pep])[,5])),2)
  POPG$meanD[i]<-round(mean(abs(read.xvg(archivos[pep])[,5])),2)
  POPG$time[i]<-length(which(abs(read.xvg(archivos[pep])[,5])<=2))*3
  POPG$ftime[i]<-((which(abs(read.xvg(archivos[pep])[,5])<=2)[1])*3)/1000
}
#Energía
setwd("/Users/Daniel/Desktop/BioInformatics_Data/")
archivos<-list.files("ENER/ENER_POPG/")
setwd("ENER/ENER_POPG/")
for (i in (1:59)){
  pep<-grep(as.vector(POPG$ID[i]),archivos)[1]
  A<-read.xvg(archivos[pep])
  POPG$ener[i]<-min(A[,2]-A[1,2])
}

#RMSD
setwd("/Users/Daniel/Desktop/BioInformatics_Data/")
archivos<-list.files("RMSD/RMSD_POPG/")
setwd("RMSD/RMSD_POPG")
for (i in c(1:60)){
  pep<-grep(as.vector(POPG$ID)[i],archivos)[1]
  POPG$RMSD[i]<-round(mean(read.xvg(archivos[pep])[,2]),2)
}

#Giro
setwd("/Users/Daniel/Desktop/BioInformatics_Data/")
archivos<-list.files("GYR/GYR_POPG/")
setwd("GYR/GYR_POPG/")
for (i in (1:60)){
  pep<-grep(as.vector(POPG$ID)[i],archivos)[1]
  POPG$GIRO[i]<-round(mean(read.xvg(archivos[pep])[,2]),2)
}

###################
# GENERACIÓN DE ARCHIVOS Y CÁLCULO DEL APL
###################

#APL
setwd("/Users/Daniel/Desktop/BioInformatics_Data/")
archivos<-list.files("MD/MD_POPG/")
setwd("MD/MD_POPG/")
for (i in c(1:60)){
  pep<-grep(as.vector(POPG$ID)[i],archivos)[1]
  pdb<-read.pdb(archivos[pep])
  write.pdb(pdb,paste("VTMC-",POPG$ID[i],".pdb",sep=""))
  sink(paste("VTMC-",POPG$ID[i],sep=""))
  writeLines("[INPUT]")
  writeLines(paste("INPUTPDB= ","L-VTMC-",POPG$ID[i],".pdb-A",sep=""))
  writeLines("\n")
  writeLines("[OUTPUT]")
  writeLines(paste("OUTNAME  = ","L-VTMC-",POPG$ID[i],".pdb-A",sep=""))
  writeLines("OUTPUTLEV = 0")
  writeLines("\n")
  writeLines("[PROTEIN]")
  SEQ<-vector()
  SEQ[2:(length(s2c(as.vector(POPG$SEQUENCE[i])))+1)]<-toupper(aaa(s2c(as.vector(POPG$SEQUENCE[i]))))
  SEQ[1]<-"SELECT  = "
  writeLines(SEQ,sep=" ")
  writeLines("\n
HYDROGEN = INCLUDE
\n
[LIPID]
SELECT = POP
HEAD_ATOM = P8
HYDROGEN = INCLUDE
\n
[ATOM_MASS]
H =  1.008
C =  12.01
O =  16.00
N =  14.01
P =  30.97
S =  32.07
\n
[VDW_RADIUS]
H =  1.2
C =  1.7
O =  1.52
N =  1.55
P =  1.8
S =  1.8
\n
[ANALYSIS]")
  writeLines(paste("CELL_SIZE =",as.character(pdb$cryst1$abc[1]),as.character(pdb$cryst1$abc[2]),sep=" "))
  writeLines("CELL_CENTER =    0.0   0.0
RADIUS      =  0.1
FACTOR      =  100
ISEED       =  3141592
RANDOMIZE     =  YES
DECIMAL_PLACE =  12")
  sink()
}

#############
# AGREGAR LAS LINEAS 73-76
# USAR ARREGLARPDB.PY
# CORRER VTMC
#############

# LECTURA DE ARCHIVOS DE SALIDA DEL APL
#############
setwd("/Users/Daniel/Desktop/BioInformatics_Data/")
archivos<-list.files("MD/MD_POPG/",pattern=".out")
setwd("MD/MD_POPG/")
for (i in c(1:60)){
  peptido<-grep(as.vector(POPG$ID[i]),archivos)[1]
  archivo<-readLines(archivos[peptido])
  linea<-grep("non-boundary",archivo)
  POPG$bAPL[i]<-as.numeric(scan(archivos[peptido],skip=linea-5,what="n",quiet=TRUE)[21])
  POPG$nLb[i]<-as.numeric(scan(archivos[peptido],skip=linea-5,what="n",quiet=TRUE)[23])
  POPG$nbAPLP[i]<-as.numeric(scan(archivos[peptido],skip=linea-5,what="n",quiet=TRUE)[30])
}

###########
# ARCHIVOS DE CONFIGURACIÓN GRIDMAT
###########
#Thickness
setwd("/Users/Daniel/Desktop/BioInformatics_Data/")
archivos<-list.files("MD/MD_POPG/",pattern=".pdb-A")
setwd("MD/MD_POPG/")
for (i in 1: length(POPG$ID)){
  pep<-grep(as.vector(POPG$ID)[i],archivos)[1]
  sink(paste("GRIDMAT-",POPG$ID[i],sep=""))
  writeLines(paste("coord_file", archivos[pep],sep=" ")) 
  writeLines("file_type pdb
num_frames 1 
num_lipid_types 1
resname1 POP 
atomname1 P8 
solvent SOL 
ions NA+,CL-
box_size solvent 
grid 100 
conserve_ratio no
protein yes")
  writeLines(paste("output_prefix", archivos[pep],sep=" ")) 
  writeLines("output_format column 
thickness yes 
area no")
  sink()
}

#############
# LECTURA DE RESULTADOS
#############
archivos<-list.files("./",pattern="bottom")
for (i in c(1:60)){
  pep<-grep(as.vector(POPG$ID[i]),archivos)[1]
  POPG$maxThick[i]<-round(max(read.csv(archivos[pep])[,1]),2)
  POPG$meanThick[i]<-round(mean(read.csv(archivos[pep])[,1]),2)
}
