package PTE;

# Copyright 2012 Thomas H. Schmidt & Christian Kandt
#
# This file is part of LAMBADA.
#
# LAMBADA is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# LAMBADA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with LAMBADA; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

use strict;

require Exporter;

my @ISA = qw(Exporter);
my @EXPORT = qw(getRadius);

my @atomicNum2Element = qw(0 'H' 'He');


my %element2atomicNum = ('H' => 1,
                         'He' => 2);

my @elements = (
    {
        name      => 'Nothing',
        symbol    => 'No',
        number    => 0,
        weight    => 0,
        covRadius => 0, 
        vdwRadius => 0
    },
    {
        name      => 'Hydrogen',
        symbol    => 'H',
        number    => 1,
        weight    => 1.007947,
        covRadius => 0.031, 
        vdwRadius => 0.11
    },
    {
        name      => 'Helium',
        symbol    => 'He',
        number    => 2,
        weight    => 4.0026022,
        covRadius => 0.028, 
        vdwRadius => 0.14
    }
);


my %elementsBySymbol;
foreach (@elements) {
    $elementsBySymbol{$$_{symbol}} = $_;
}



sub addElement {
    push(@elements, $_[0]);
}



sub getSymbol {
    return $elements[$_[0]]{symbol};
}



sub getName {
    return $elements[$_[0]]{name};
}



sub getCovRadius {
    return $elements[$_[0]]{covRadius};
}



sub getVdwRadius {
    return $elements[$_[0]]{vdwRadius};
}



sub getRadius {
    my %atomRadius = ('H'  => 0.110,
                      'He' => 0.140,
                      'Li' => 0.182,
                      'Be' => 0.153,
                      'B'  => 0.192,
                      'C'  => 0.170,
                      'N'  => 0.155,
                      'O'  => 0.152,
                      'F'  => 0.147,
                      'Ne' => 0.154,
                      'Na' => 0.227,
                      'Mg' => 0.173,
                      'Al' => 0.184,
                      'Si' => 0.210,
                      'P'  => 0.180,
                      'S'  => 0.180,
                      'Cl' => 0.175,
                      'Ar' => 0.188,
                      'K'  => 0.275,
                      'Ca' => 0.231,
                      'Ni' => 0.163,
                      'Cu' => 0.140,
                      'Zn' => 0.139,
                      'Ga' => 0.187,
                      'Ge' => 0.211,
                      'As' => 0.185,
                      'Se' => 0.190,
                      'Br' => 0.185,
                      'Kr' => 0.202,
                      'Rb' => 0.303,
                      'Sr' => 0.249,
                      'Pd' => 0.163,
                      'Ag' => 0.172,
                      'Cd' => 0.158,
                      'In' => 0.193,
                      'Sn' => 0.217,
                      'Sb' => 0.206,
                      'Te' => 0.206,
                      'I'  => 0.198,
                      'Xe' => 0.216,
                      'Cs' => 0.343,
                      'Ba' => 0.268,
                      'Pt' => 0.175,
                      'Au' => 0.166,
                      'Hg' => 0.155,
                      'Tl' => 0.196,
                      'Pb' => 0.202,
                      'Bi' => 0.207,
                      'Po' => 0.197,
                      'At' => 0.202,
                      'Rn' => 0.220,
                      'Fr' => 0.348,
                      'Ra' => 0.283,
                      'U'  => 0.186,
                      'CGA'=> 0.400);

    return $atomRadius{$_[0]} if $atomRadius{$_[0]};
    return 0;

# Atomic van der Waals radii in nm.
#    1) R. Scott Rowland, Robin Taylor: Intermolecular Nonbonded Contact
#       Distances in Organic Crystal Structures: Comparison with Distances
#       Expected from van der Waals Radii. In: J. Phys. Chem. 1996, 100,
#       S. 7384–7391, doi:10.1021/jp953141+.
#    2) A. Bondi: van der Waals Volumes and Radii. In: J. Phys. Chem. 1964, 68,
#       S. 441-451, doi:10.1021/j100785a001.
#    3) Manjeera Mantina, Adam C. Chamberlin, Rosendo Valero, Christopher J.
#       Cramer, Donald G. Truhlar: Consistent van der Waals Radii for the Whole
#       Main Group. In: J. Phys. Chem. A. 2009, 113, S. 5806–5812,
#       doi:10.1021/jp8111556.
}

1;