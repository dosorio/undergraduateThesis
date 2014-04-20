#!/bin/sh

#  InterPEPT.sh
#  SCRIPT PARA EVALUAR LA INTERACCIÓN DE PÉPTIDOS CON BICAPAS LIPÍDICAS
#  DANIEL CAMILO OSORIO
#  UNIVERSIDAD INDUSTRIAL DE SANTANDER
#  d.osorio@me.com



# PATH
#export GMX_MAXCONSTRWARN=-1
export PATH=$PATH:/usr/bin:/usr/local/cuda/bin:/usr/local/gromacs/bin
export DSSP=/home/dosorio/DSSP
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64:/usr/local/gromacs/lib

# CREANDO CARPETA DE RESULTADOS
mkdir PEP-MEM ; cd PEP-MEM

# DEFINIENDO PÉPTIDOS A SIMULAR
for pept in 2JQ0 2M9I GIB10 GIB20 GIB30 1KUW 1OT0 2K10 2RSH GIB1 GIB11 GIB21 1S6W 2K60 2KAM 3Q8J GIB2 GIB12 GIB22 1T51 2KHF 2LL1 4B19 GIB3 GIB13 GIB23 1WO0 2LNF 4B2U GIB4 GIB14 GIB24 2MAG 1X22 2LO7 4BMF GIB5 GIB15 GIB25 2RLH 1Z64 2LQ0 GIB6 GIB16 GIB26 2AMN 2LQ1 GIB7 GIB17 GIB27 2AP7 2LQA GIB8 GIB18 GIB28 2B68 2M0D GIB9 GIB19 GIB29
do

# CREANDO LA CARPETA DEL PÉPTIDO
mkdir $pept ; cd $pept

for memb in POPG POPE POPC
do

# CREANDO CARPETA EN RESULTADOS
mkdir $memb ; cd $memb

# INCLUYE LAS CONFIGURACIONES Y ESTRUCTURAS
cp ~/PDB/$pept.pdb ./
cp ~/MEM/$memb* ./
cp -r ~/FF/* ./

if [ $memb = POPG ];
then
cp ~/MEM/dPOPG.itp ./
fi

#CREAR LA TOPOLOGIA
pdb2gmx_mpi -f $pept.pdb -o $pept.gro -ignh -water spc<<EOF
14
EOF

#ELIMINAR LA PERIODICIDAD DE LA MEMBRANA
#editconf -f $membrana.pdb -o $membrana.gro -pbc

# CREANDO EL ARCHIVO MINIM.MDP
cat > minim.mdp << EOF
integrator	= steep
emtol		= 1000.0
emstep      = 0.001
nsteps		= 500000
nstlist		= 1
ns_type		= grid
rlist		= 1.2
coulombtype	= PME
rcoulomb	= 1.2
rvdw		= 1.2
pbc         = xyz
EOF

# CONFIGURAR EL T=0 DE LA MEMBRANA
grompp_mpi -f minim.mdp -c $memb.gro -p $memb.top -o em.tpr

# CREAR LA CAJA DE SIMULACIÃ“N DE LA MEMBRANA
trjconv_mpi -s em.tpr -f $memb.gro -o n-$memb.gro -pbc mol -ur compact <<EOF
2
EOF

# AUMENTANDO TAMAÑO EN EL EJE Z
if [ $memb = POPC ];
then
editconf_mpi -f n-$memb.gro -o w-$memb.gro -c -box 6.48650   6.48650    9.00
else
if [ $memb = POPE ];
then
editconf_mpi -f n-$memb.gro -o w-$memb.gro -c -box 9.57070   9.48610    9.00
else
if [ $memb = POPG ];
then
editconf_mpi -f n-$memb.gro -o w-$memb.gro -c -box 6.63780   6.62730    9.00
fi
fi
fi

# POSICIONANDO EL PEPTIDO
if [ $memb = POPC ];
then
editconf_mpi -f $pept.gro -o nb-$pept.gro -box 6.48650   6.48650   9.00000 -center 3.24325  3.24325  0.5
else
if [ $memb = POPE ];
then
editconf_mpi -f $pept.gro -o nb-$pept.gro -box 9.57070   9.48610   9.00000 -center 4.78535  4.74305  0.5
else
if [ $memb = POPG ];
then
editconf_mpi -f $pept.gro -o nb-$pept.gro -box 6.63780   6.62730   9.00000 -center 3.31800  3.31300  0.5
fi
fi
fi

# CREANDO EL ARCHIVO DE UNIÓN DE TOPOLOGÍAS
cp ~/CONFIG/UnirGRO.py ./

#CREAR LA CAJA CONJUNTA DE SIMULACIÃ“N
cat nb-$pept.gro w-$memb.gro > system.gro
python UnirGRO.py
rm system.gro
mv system2.gro $pept-$memb.gro

# UNIENDO TOPOLOGÍAS
cp ~/CONFIG/$memb.py ./
python $memb.py

# RENOMBRAMIENTO DE ARCHIVOS
rm topol.top
mv topol2.top topol.top

# CREACIÓN DEL ARCHIVO VDWRADII PARA SOLVATACIÓN
cat > vdwradii.dat <<EOF
; Very approximate VanderWaals radii
; only used for drawing atoms as balls or for calculating atomic overlap.
; longest matches are used
; '???' or '*' matches any residue name
; 'AAA' matches any protein residue name
???  C     0.375
???  F     0.12
???  H     0.04
???  N     0.110
???  O     0.105
???  S     0.16
; Water charge sites
SOL  MW    0
SOL  LP    0
; Masses for vsite construction
GLY  MN1   0
GLY  MN2   0
ALA  MCB1  0
ALA  MCB2  0
VAL  MCG1  0
VAL  MCG2  0
ILE  MCG1  0
ILE  MCG2  0
ILE  MCD1  0
ILE  MCD2  0
LEU  MCD1  0
LEU  MCD2  0
MET  MCE1  0
MET  MCE2  0
TRP  MTRP1 0
TRP  MTRP2 0
THR  MCG1  0
THR  MCG2  0
LYSH MNZ1  0
LYSH MNZ2  0
EOF

#AÃ‘ADIR MOLECULAS DE AGUA
genbox_mpi -cp $pept-$memb.gro -cs spc216.gro -o $pept-$memb-solv.gro -p topol.top

#ELIMINAR EL ;VDW
rm vdwradii.dat

# CREANDO EL ARCHIVO ions.mdp
cat > ions.mdp <<EOF
integrator	= steep
emtol		= 1000.0
emstep 	        = 0.01
nsteps		= 50000
nstlist		= 1
ns_type		= grid
rlist		= 1.2
coulombtype	= PME
rcoulomb	= 1.2
rvdw		= 1.2
pbc		= xyz
EOF

#CONFIGURAR LA ADICIÃ“N DE IONES
grompp_mpi -f ions.mdp -c $pept-$memb-solv.gro -p topol.top -o ions.tpr


#AÃ‘ADIR IONES HASTA NEUTRALIZAR LA CARGA
genion_mpi -s ions.tpr -o i-$pept-$memb.gro -p topol.top -pname NA -nname CL -conc 0.9 -neutral -rmin 0.1<<EOF
15
EOF

#CONFIGURACIÃ“N DE LA MINIMIZACIÃ“N DE ENERGÃA
grompp_mpi -f minim.mdp -c i-$pept-$memb.gro -p topol.top -o $pept-$memb.tpr

#MINIMIZACIÃ“N HASTA 1000KJ/MOL
mpirun -np 8 mdrun_mpi -deffnm $pept-$memb

#HACER GRUPOS DE SIMULACIÃ“N
make_ndx_mpi -f $pept-$memb.gro -o index-$pept-$memb.ndx <<EOF
17 | 19
1 | 13
q
EOF

# SIMULACIONES A DIFERENTES TEMPERATURAS
for temp in 293 310 323
do
if [ $memb = POPG ];
then
cat > nvt.mdp <<EOF
title					= NVT
integrator				= md
nsteps					= 50000
dt	   			 	= 0.002
nstxout					= 100
nstvout					= 100
nstenergy				= 100
nstlog					= 100
continuation				= no
constraint_algorithm 	= lincs
constraints				= all-bonds
lincs_iter				= 1
lincs_order				= 4
ns_type					= grid
nstlist					= 5
rlist					= 1.2
rcoulomb				= 1.2
rvdw					= 1.2
coulombtype				= PME
pme_order				= 4
fourierspacing			= 0.16
tcoupl					= V-rescale
tc-grps					= Protein DPOPG SOL_Ion
tau_t					= 0.1	0.1	0.1
ref_t					= $temp $temp $temp
pcoupl					= no
pbc		    			= xyz
DispCorr				= EnerPres
gen_vel					= yes
gen_temp				= $temp
gen_seed				= -1
nstcomm					= 1
comm-mode				= Linear
comm-grps				= Protein_DPOPG SOL_Ion
cutoff-scheme           = Verlet
EOF
else
cat > nvt.mdp <<EOF
title					= NVT
integrator				= md
nsteps					= 50000
dt		   			= 0.002
nstxout					= 100
nstvout					= 100
nstenergy				= 100
nstlog					= 100
continuation			= no
constraint_algorithm 	= lincs
constraints				= all-bonds
lincs_iter				= 1
lincs_order				= 4
ns_type					= grid
nstlist					= 5
rlist					= 1.2
rcoulomb				= 1.2
rvdw					= 1.2
coulombtype				= PME
pme_order				= 4
fourierspacing			= 0.16
tcoupl					= V-rescale
tc-grps					= Protein $memb SOL_Ion
tau_t					= 0.1	0.1	0.1
ref_t					= $temp $temp $temp
pcoupl					= no
pbc		    			= xyz
DispCorr				= EnerPres
gen_vel					= yes
gen_temp				= $temp
gen_seed				= -1
nstcomm					= 1
comm-mode				= Linear
comm-grps				= Protein_$memb SOL_Ion
cutoff-scheme           = Verlet
EOF
fi

#ACOPLAMIENTO NVT EN X TEMPERATURA
grompp_mpi -f nvt.mdp -c $pept-$memb.gro -p topol.top -n index-$pept-$memb.ndx -o nvt-$pept-$memb-$temp.tpr
mpirun -np 8 mdrun_mpi -deffnm nvt-$pept-$memb-$temp 

# CREANDO EL ARCHIVO npt.mdp
if [ $memb = POPG ]
then
cat > npt.mdp <<EOF
integrator          = md
nsteps              = 50000
dt                  = 0.002
nstxout             = 100
nstvout             = 100
nstenergy           = 100
nstlog              = 100
continuation        = yes
constraint_algorithm = lincs
constraints         = all-bonds
lincs_iter          = 1
lincs_order         = 4
ns_type             = grid
nstlist             = 5
rlist               = 1.2
rcoulomb            = 1.2
rvdw                = 1.2
coulombtype         = PME
pme_order           = 4
fourierspacing      = 0.16
tcoupl              = Nose-Hoover
tc-grps             = Protein DPOPG	SOL_Ion
tau_t               = 0.5	0.5	0.5
ref_t               = $temp $temp $temp
pcoupl              = Parrinello-Rahman
pcoupltype          = semiisotropic
tau_p               = 5.0
ref_p               = 1.0	1.0
compressibility     = 4.5e-5	4.5e-5
pbc                 = xyz
DispCorr            = EnerPres
gen_vel             = no
cutoff-scheme       = Verlet
EOF
else
cat > npt.mdp <<EOF
integrator          = md
nsteps              = 50000
dt                  = 0.002
nstxout             = 100
nstvout             = 100
nstenergy           = 100
nstlog              = 100
continuation        = yes
constraint_algorithm = lincs
constraints         = all-bonds
lincs_iter          = 1
lincs_order         = 4
ns_type             = grid
nstlist             = 5
rlist               = 1.2
rcoulomb            = 1.2
rvdw                = 1.2
coulombtype         = PME
pme_order           = 4
fourierspacing      = 0.16
tcoupl              = Nose-Hoover
tc-grps             = Protein $memb	SOL_Ion
tau_t               = 0.5	0.5	0.5
ref_t               = $temp $temp $temp
pcoupl              = Parrinello-Rahman
pcoupltype          = semiisotropic
tau_p               = 5.0
ref_p               = 1.0	1.0
compressibility     = 4.5e-5	4.5e-5
pbc                 = xyz
DispCorr            = EnerPres
gen_vel             = no
cutoff-scheme       = Verlet
EOF
fi

#ACOPLAMIENTO NPT EN X TEMPERATURA
grompp_mpi -f npt.mdp -c nvt-$pept-$memb-$temp.gro -t nvt-$pept-$memb-$temp.cpt -p topol.top -n index-$pept-$memb.ndx -o npt-$pept-$memb-$temp.tpr
mpirun -np 8 mdrun_mpi -deffnm npt-$pept-$memb-$temp 

if [ $memb = POPG ]
then
# CREANDO EL ARCHIVO md.mdp
cat > md.mdp <<EOF
integrator              = md
nsteps                  = 16666667
dt                      = 0.003
nstxout                 = 1000
nstvout                 = 1000
nstxtcout               = 1000
nstenergy               = 1000
nstlog                  = 1000
continuation            = yes
constraint_algorithm    = lincs
constraints             = all-bonds
lincs_iter              = 1
lincs_order             = 4
ns_type                 = grid
nstlist                 = 5
rlist                   = 1.2
rcoulomb                = 1.2
rvdw                    = 1.2
coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16
tcoupl                  = Nose-Hoover
tc-grps                 = Protein DPOPG	SOL_Ion
tau_t                   = 0.5	0.5	0.5
ref_t                   = $temp $temp $temp
pcoupl                  = Parrinello-Rahman
pcoupltype              = semiisotropic
tau_p                   = 2.0
ref_p                   = 1.0	1.0
compressibility         = 4.5e-5	4.5e-5
pbc                     = xyz
DispCorr                = EnerPres
gen_vel                 = no
cutoff-scheme           = Verlet
energygrps		= Protein DPOPG
;energygrps-table        = Protein DPOPG
EOF
else
cat > md.mdp <<EOF
integrator              = md
nsteps                  = 16666667
dt                      = 0.003
nstxout                 = 1000
nstvout                 = 1000
nstxtcout               = 1000
nstenergy               = 1000
nstlog                  = 1000
continuation            = yes
constraint_algorithm    = lincs
constraints             = all-bonds
lincs_iter              = 1
lincs_order             = 4
ns_type                 = grid
nstlist                 = 5
rlist                   = 1.2
rcoulomb                = 1.2
rvdw                    = 1.2
coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16
tcoupl                  = Nose-Hoover
tc-grps                 = Protein $memb	SOL_Ion
tau_t                   = 0.5	0.5	0.5
ref_t                   = $temp $temp $temp
pcoupl                  = Parrinello-Rahman
pcoupltype              = semiisotropic
tau_p                   = 2.0
ref_p                   = 1.0	1.0
compressibility         = 4.5e-5	4.5e-5
pbc                     = xyz
DispCorr                = EnerPres
gen_vel                 = no
cutoff-scheme           = Verlet
energygrps              = Protein $memb
;energygrps-table        = Protein $memb
EOF
fi


#DINÃMICA MOLECULAR
grompp_mpi -f md.mdp -c npt-$pept-$memb-$temp.gro -t npt-$pept-$memb-$temp.cpt -p topol.top -n index-$pept-$memb.ndx -o md-$pept-$memb-$temp.tpr
mpirun -np 8 mdrun_mpi -deffnm md-$pept-$memb-$temp -v

# GRAFICO DE MINIMIZACIÃ“N DE ENERGÃA
g_energy_mpi -f $pept-$memb.edr -o R_epot-$pept-$memb.xvg<<EOF
11
0
EOF

# GRÃFICO DE TEMPERATURA
g_energy_mpi -f nvt-$pept-$memb-$temp.edr -o R_temp-$pept-$temp-$memb.xvg<<EOF
15
0
EOF

# GRÃFICO DE PRESIÃ“N
g_energy_mpi -f npt-$pept-$memb-$temp.edr -o R_pres-$pept-$temp-$memb.xvg<<EOF
16
0
EOF

# GRÃFICO DE DENSIDAD
g_energy_mpi -f npt-$pept-$memb-$temp.edr -o R_den-$pept-$temp-$memb.xvg<<EOF
22
0
EOF

# CONVIRTIENDO LAS TRAYECTORIAS
trjconv_mpi -s md-$pept-$memb-$temp.tpr -f md-$pept-$memb-$temp.xtc -o nopbc-$pept-$memb-$temp.xtc -pbc mol -ur compact <<EOF
0
EOF

# RMSD
g_rms_mpi -s $pept-$memb.tpr -f nopbc-$pept-$memb-$temp.xtc -o R_rmsd-$pept-$memb-$temp.xvg -tu ns <<EOF
4
4
EOF

# GIRO
g_gyrate_mpi -s md-$pept-$memb-$temp.tpr -f nopbc-$pept-$memb-$temp.xtc -o R_giro-$pept-$memb-$temp.xvg<<EOF
1
EOF

# CONVIRTIENDO TRAYECTORIAS
trjconv_mpi -f md-$pept-$memb-$temp.xtc -o R_md-$pept-$memb-$temp-1.pdb -s md-$pept-$memb-$temp.tpr -n index-$pept-$memb.ndx -pbc mol -b 4991 -e 5009 <<EOF
0
EOF
trjconv_mpi -f md-$pept-$memb-$temp.xtc -o R_md-$pept-$memb-$temp-2.pdb -s md-$pept-$memb-$temp.tpr -n index-$pept-$memb.ndx -pbc mol -b 9991 -e 10009 <<EOF
0
EOF
trjconv_mpi -f md-$pept-$memb-$temp.xtc -o R_md-$pept-$memb-$temp-3.pdb -s md-$pept-$memb-$temp.tpr -n index-$pept-$memb.ndx -pbc mol -b 14991 -e 15009 <<EOF
0
EOF
trjconv_mpi -f md-$pept-$memb-$temp.xtc -o R_md-$pept-$memb-$temp-4.pdb -s md-$pept-$memb-$temp.tpr -n index-$pept-$memb.ndx -pbc mol -b 19991 -e 20009 <<EOF
0
EOF
trjconv_mpi -f md-$pept-$memb-$temp.xtc -o R_md-$pept-$memb-$temp-5.pdb -s md-$pept-$memb-$temp.tpr -n index-$pept-$memb.ndx -pbc mol -b 24991 -e 25009 <<EOF
0
EOF
trjconv_mpi -f md-$pept-$memb-$temp.xtc -o R_md-$pept-$memb-$temp-6.pdb -s md-$pept-$memb-$temp.tpr -n index-$pept-$memb.ndx -pbc mol -b 29991 -e 30009 <<EOF
0
EOF
trjconv_mpi -f md-$pept-$memb-$temp.xtc -o R_md-$pept-$memb-$temp-7.pdb -s md-$pept-$memb-$temp.tpr -n index-$pept-$memb.ndx -pbc mol -b 34991 -e 35009 <<EOF
0
EOF
trjconv_mpi -f md-$pept-$memb-$temp.xtc -o R_md-$pept-$memb-$temp-8.pdb -s md-$pept-$memb-$temp.tpr -n index-$pept-$memb.ndx -pbc mol -b 39991 -e 40009 <<EOF
0
EOF
trjconv_mpi -f md-$pept-$memb-$temp.xtc -o R_md-$pept-$memb-$temp-9.pdb -s md-$pept-$memb-$temp.tpr -n index-$pept-$memb.ndx -pbc mol -b 44991 -e 45009 <<EOF
0
EOF
trjconv_mpi -f md-$pept-$memb-$temp.xtc -o R_md-$pept-$memb-$temp-10.pdb -s md-$pept-$memb-$temp.tpr -n index-$pept-$memb.ndx -pbc mol -b 49982 -e 50000 <<EOF
0
EOF


# E3D - DSSP
#do_dssp_mpi -f nopbc-$pept-$memb-$temp.xtc -s md-$pept-$memb-$temp.tpr -n index-$pept-$memb.ndx -o R_E3D-$pept-$memb-$temp.xpm -tu ns <<EOF
#1
#EOF

# APL
g_energy_mpi -f md-$pept-$memb-$temp.edr -o R_apl-$pept-$memb-$temp.xvg <<EOF
18 19
0
EOF

# AGUA
trjconv_mpi -f md-$pept-$memb-$temp.xtc -o H2O-$pept-$memb-$temp.pdb -s md-$pept-$memb-$temp.tpr -n index-$pept-$memb.ndx -pbc mol << EOF
16
EOF
awk '{print $8}' H2O-$pept-$memb-$temp.pdb > R_H2O-$pept-$memb-$temp.txt

# DISTANCIA PEPT-PROTEINA
g_dist_mpi -f md-$pept-$memb-$temp.xtc -s md-$pept-$memb-$temp.tpr -n index-$pept-$memb.ndx -o R_dist-$pept-$memb-$temp.xvg << EOF
1
13
EOF

# ENERGÍA INTERACCIÓN PEPTIDO PROTEINA
g_energy_mpi -f md-$pept-$memb-$temp.edr -s md-$pept-$memb-$temp.tpr -o R_ener-$pept-$memb-$temp.xvg <<EOF
6 9 10 51 53 0
EOF

done
cp -r R_* ~/R_INTER/
rm -r *
cd ..
done
cd ..
done
exit

