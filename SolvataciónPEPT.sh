#!/bin/sh

#  SOLV.sh
#  SCRIPT PARA EVALUAR LA ESTABILIDAD DE LOS PÉPTIDOS CATIÓNICOS ANTIMICROBIANOS EN AGUA
#  DANIEL CAMILO OSORIO
#  UNIVERSIDAD INDUSTRIAL DE SANTANDER
#  d.osorio@me.com

# PATH
export PATH=$PATH:/usr/bin:/usr/local/cuda/bin:/usr/local/gromacs/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64:/usr/local/gromacs/lib

# DEFINIENDO PÉPTIDOS A SIMULAR
for pept in 1KUW  1WO0  2AP7  2K6O  2LNF  2LQA  2RLH  4B2U  1OT0  1X22  2B68  2KAM  2LO7  2M0D  2RSH  4BMF 1S6W  1Z64  2JQ0  2KHF  2LQ0  2M9I  3Q8J 1T51  2AMN  2K10  2LL1  2LQ1  2MAG  4B19
do

# CREANDO CARPETAS DE SIMULACIÓN
mkdir ~/SOLV
mkdir ~/SOLV/$pept
cp ~/E/$pept.pdb ~/SOLV/$pept/
cd ~/SOLV/$pept/

# GENERA EL ARCHIVO DE ADICIÓN DE IONES
cat > ions.mdp << EOF
integrator      = steep
emtol		    = 1000.0
emstep          = 0.01
nsteps          = 50000
nstlist         = 1
ns_type         = grid
rlist		    = 1.0
coulombtype     = PME
rcoulomb        = 1.0
rvdw		    = 1.0
pbc             = xyz
EOF

# GENERA EL ARCHIVO DE MINIMIZACIÓN DE ENERGÍA
cat > minim.mdp << EOF
integrator	= steep
emtol		= 1000.0
emstep      = 0.01
nsteps		= 50000
nstlist		= 1
ns_type		= grid
rlist		= 1.0
coulombtype	= PME
rcoulomb	= 1.0
rvdw		= 1.0
pbc         = xyz
EOF

# CONVERSIÓN DEL PDB AL GRO
pdb2gmx_mpi -f $pept.pdb -o $pept.gro -ignh <<EOF
13	# Utiliza el Campo de Fuerza Gromos53a6
1	# SPC Modelo de Agua
EOF


# DEFINIR EL TAMAÑO DE LA CAJA DE SIMULACIÓN
editconf_mpi -f $pept.gro -o $pept-NB.gro -c -d 1.0 -bt cubic

# AÑADIR MOLECULAS DE AGUA
genbox_mpi -cp $pept-NB.gro -cs spc216.gro -o $pept-SOLV.gro -p topol.top

# CONFIGURAR LA ADICIÓN DE IONES
grompp_mpi -f ions.mdp -c $pept-SOLV.gro -p topol.top -o ions.tpr

# AÑADIR IONES
genion_mpi -s ions.tpr -o $pept-ION.gro -p topol.top -pname NA -nname CL -neutral <<EOF
13 # Añade los iones necesarios al solvente para estabilizar el sistema
EOF

# MINIMIZACIÓN DE ENERGÍA
if [ $pept = 1KUW ] || [ $pept = 2M0D ] || [ $pept = 4B2U ] || [ $pept = 4BMF ];
then
grompp_mpi -f minim.mdp -c $pept-SOLV.gro -p topol.top -o em-$pept.tpr
mpirun -np 8 mdrun_mpi -v -deffnm em-$pept
else
grompp -f minim.mdp -c $pept-ION.gro -p topol.top -o em-$pept.tpr
mpirun -np 8 mdrun_mpi -v -deffnm em-$pept
fi

# GRAFICO DE MINIMIZACIÓN DE ENERGÍA
g_energy_mpi -f em-$pept.edr -o epot-$pept.xvg<<EOF
10 0
EOF

for temp in 293 310 323
do

cat > nvt.mdp << EOF
integrator              = md
nsteps                  = 50000
dt                      = 0.002
nstxout                 = 100
nstvout                 = 100
nstenergy               = 100
nstlog                  = 100
continuation            = no
constraint_algorithm    = lincs
constraints             = all-bonds
lincs_iter              = 1
lincs_order             = 4
ns_type                 = grid
nstlist                 = 5
rlist                   = 1.0
rcoulomb                = 1.0
rvdw                    = 1.0
coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16
tcoupl                  = V-rescale
tc-grps                 = Protein Non-Protein
tau_t                   = 0.1	0.1
ref_t                   = $temp $temp
pcoupl                  = no
pbc                     = xyz
DispCorr                = EnerPres
gen_vel                 = yes
gen_temp                = 293
gen_seed                = -1
cutoff-scheme           = Verlet
EOF

cat > npt.mdp << EOF
integrator              = md
nsteps                  = 50000
dt                      = 0.002
nstxout                 = 100
nstvout                 = 100
nstenergy               = 100
nstlog                  = 100
continuation            = yes
constraint_algorithm    = lincs
constraints             = all-bonds
lincs_iter              = 1
lincs_order             = 4
ns_type                 = grid
nstlist                 = 5
rlist                   = 1.0
rcoulomb                = 1.0
rvdw                    = 1.0
coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16
tcoupl                  = V-rescale
tc-grps                 = Protein Non-Protein
tau_t                   = 0.1	0.1
ref_t                   = $temp $temp
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5
refcoord_scaling        = com
pbc                     = xyz
DispCorr                = EnerPres
gen_vel                 = no
cutoff-scheme           = Verlet
EOF

cat > md.mdp << EOF
integrator              = md
nsteps                  = 1000000
dt                      = 0.002
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
rlist                   = 1.0
rcoulomb                = 1.0
rvdw                    = 1.0
coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16
tcoupl                  = V-rescale
tc-grps                 = Protein Non-Protein
tau_t                   = 0.1	0.1
ref_t                   = $temp $temp
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5
pbc                     = xyz
DispCorr                = EnerPres
gen_vel                 = no
cutoff-scheme           = Verlet
EOF

# ACOPLAMIENTO NVT
grompp_mpi -f nvt.mdp -c em-$pept.gro -p topol.top -o nvt-$pept-$temp.tpr
mpirun -np 8 mdrun_mpi -v -deffnm nvt-$pept-$temp

# GRÁFICO DE TEMPERATURA
g_energy_mpi -f nvt-$pept-$temp.tpr -o $temp-$pept-temp.xvg<<EOF
15 0
EOF

# ACOPLAMIENTO NPT
grompp_mpi -f npt.mdp -c nvt-$pept-$temp.gro -t nvt-$pept-$temp.cpt -p topol.top -o npt-$pept-$temp.tpr
mdrun -v -deffnm npt-$pept-$temp

# GRÁFICO DE PRESIÓN
g_energy_mpi -f npt-$pept-$temp.edr -o pres-$pept-$temp.xvg<<EOF
16 0
EOF

# GRÁFICO DE DENSIDAD
g_energy_mpi -f npt-$pept-$temp.edr -o den-$pept-$temp.xvg<<EOF
22 0
EOF

# DINÁMICA MOLECULAR
grompp_mpi -f md.mdp -c npt-$pept-$temp.gro -t npt-$pept-$temp.cpt -p topol.top -o md-$pept-$temp.tpr
mpirun -np 8 mdrun_mpi -v -deffnm md-$pept-$temp

# CONVIRTIENDO LAS TRAYECTORIAS
trjconv_mpi -s md-$pept-$temp.tpr -f md-$pept-$temp.xtc -o nopbc-$pept-$temp.xtc -pbc mol -ur compact <<EOF
0
EOF

# RMSD
g_rms_mpi -s em-$pept.tpr -f nopbc-$pept-$temp.xtc -o rmsd-$pept-$temp.xvg -tu ns <<EOF
4
EOF

# GIRO
g_girate_mpi -s md-$pept-$temp.tpr -f nopbc-$pept-$temp.xtc -o giro-$pept-$temp.vxg

# ELIMINANDO ARCHIVOS INTERMEDIOS
rm *.log *.mdp* *.log.*
done
cd $SOLV
done
