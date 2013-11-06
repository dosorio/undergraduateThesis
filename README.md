cAMP-Membrana
=============

**Interacciones de cAMP's con Membranas**

Código para simulaciones en GROMACS 4.6.3 entre péptidos catiónicos antimicrobianos y bicapas lípidicas puras de POPE, POPG y POPC con la finalidad de identificar secuencias biológicamente activas
```
cAMP-Membrana/
├── ESTRUCTURAS
│   └── POPULUS  **E3D de los péptidos a simular obtenidos a través del servidor 3D-JIGSAW**
│       ├── 1KUW.pdb
│       ├── 1OT0.pdb
│       ├── 1S6W.pdb
│       ├── 1T51.pdb
│       ├── 1WO0.pdb
│       ├── 1X22.pdb
│       ├── 1Z64.pdb
│       ├── 2AMN.pdb
│       ├── 2AP7.pdb
│       ├── 2B68.pdb
│       ├── 2JQ0.pdb
│       ├── 2K10.pdb
│       ├── 2K6O.pdb
│       ├── 2KAM.pdb
│       ├── 2KHF.pdb
│       ├── 2LL1.pdb
│       ├── 2LNF.pdb
│       ├── 2LO7.pdb
│       ├── 2LQ0.pdb
│       ├── 2LQ1.pdb
│       ├── 2LQA.pdb
│       ├── 2M0D.pdb
│       ├── 2M9I.pdb
│       ├── 2MAG.pdb
│       ├── 2RLH.pdb
│       ├── 2RSH.pdb
│       ├── 3Q8J.pdb
│       ├── 4B19.pdb
│       ├── 4B2U.pdb
│       └── 4BMF.pdb
├── FF  **Campo de fuerza Gromos53a6 incluyendo los lípidos de Berger**
│   ├── aminoacids.c.tdb
│   ├── aminoacids.hdb
│   ├── aminoacids.n.tdb
│   ├── aminoacids.r2b
│   ├── aminoacids.rtp
│   ├── aminoacids.vsd
│   ├── ff_dum.itp
│   ├── ffbonded.itp
│   ├── ffnonbonded.itp
│   ├── forcefield.doc
│   ├── forcefield.itp
│   ├── ions.itp
│   ├── lipid.itp
│   ├── spc.itp
│   └── watermodels.dat
├── MEMBRANAS **Configuraciones y coordenadas de las membranas a simular**
│   ├── POPC.gro
│   ├── POPC.itp
│   ├── POPC.top
│   ├── POPE.gro
│   ├── POPE.itp
│   ├── POPE.top
│   ├── POPG.gro
│   ├── POPG.itp
│   └── POPG.top
├── R **Funciones de R para cálculo de propiedades fisicoquímicas**
│   ├── AL.R
│   ├── CN.R
│   ├── II-Guruprasad.R
│   ├── Kyte-Doolittle.R
│   └── pI.R
├── README.md
├── SolvataciónPEPT.sh
├── UnirGRO.py
└── fixTOPOL.py
``` 
##Autor##
**Daniel Camilo Osorio** |
Estudiante de Biología | GIBIM | Universidad Industrial de Santander 

**Contacto:**
daniel.osorio@correo.uis.edu.co


