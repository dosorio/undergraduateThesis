package GRO;

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

our $VERSION = 1.0;
our @ISA     = qw(Exporter);
our @EXPORT  = qw(readGro writeGro renumAtoms);


sub readGro {
    my $groFile = shift;
    my %groData;

    print "  ---------------------------------\n  Read GRO file \"$groFile\"...\r";
    open(GROFILE, "<$groFile") || die "ERROR: Cannot open GRO file \"$groFile\": $!\n";
    readHeader(\*GROFILE, \%groData);
    readCoords(\*GROFILE, \%groData);
    readFooter(\*GROFILE, \%groData);
    close(GROFILE);
    print "  Read GRO file \"$groFile\": Finished\n  ---------------------------------\n\n";

    return %groData;
}



sub renumAtoms {
    my @renumAtoms;
    my $atomId = 0;
    foreach (@{$_[0]}) {
        next unless $$_{'residue'};
        $renumAtoms[++$atomId] = $_;
    }
    return \@renumAtoms;
}



sub readHeader {
    my $fileHandle = shift;
    my $groDataRef = shift;

    $$groDataRef{'title'}  = <$fileHandle>;
    $$groDataRef{'title'}  =~ s/(^\s*)|(\s*$)//g;
    $$groDataRef{'nAtoms'} = <$fileHandle>;
    $$groDataRef{'nAtoms'} =~ s/\s//g;

    print "\n    Number of atoms: " . $$groDataRef{'nAtoms'} . "\n" if $main::verbose;
}



sub readFooter {
    my $fileHandle = shift;
    my $groDataRef = shift;

    $$groDataRef{'footline'} = <$fileHandle>;
    $$groDataRef{'footline'} =~ s/(^\s*)|(\s*$)//g;
    if ($$groDataRef{'footline'} =~ /([-+]?\d*\.?\d+([eE][-+]?\d+)?)\s+([-+]?\d*\.?\d+([eE][-+]?\d+)?)\s+([-+]?\d*\.?\d+([eE][-+]?\d+)?)/) {
        $$groDataRef{'box'}{'cooX'} = $1;
        $$groDataRef{'box'}{'cooY'} = $3;
        $$groDataRef{'box'}{'cooZ'} = $5;

        print "\n    Boxsize: x=$1, y=$3, z=$5\n" if $main::verbose;
    }
}



sub readCoords {
    my $fileHandle   = shift;
    my $groDataRef   = shift;
    my $atomId       = 0; # Unique atom ID.
    my $residue      = 0; # Unique residue ID (initially 0 -> starts with 1).
    my $lastResId    = 0;

    while (<$fileHandle>) {
        chomp($_);
        unless ($_ =~ /^\s*$/) {
            $$groDataRef{'atoms'}[++$atomId] = getAtomdata($_);
            print "    Read coordinates:  $atomId\r" if $main::verbose;
            $residue++ unless $lastResId == $$groDataRef{'atoms'}[$atomId]{'resId'};
            $$groDataRef{'atoms'}[$atomId]{'residue'} = $residue;
            $lastResId = $$groDataRef{'atoms'}[$atomId]{'resId'};
        }
        return 1 if $atomId == $$groDataRef{'nAtoms'};
    }
    return 0; # If file ends before all atoms were count.
}



sub getAtomdata {
    my $atomStr = shift;
    my $strLen  = length($atomStr);
    my %atomData;

    $atomData{'resId'}    = checkSubstr($atomStr, $strLen,  0, 5);
    $atomData{'resName'}  = checkSubstr($atomStr, $strLen,  5, 5);
    $atomData{'atomName'} = checkSubstr($atomStr, $strLen, 10, 5);
    $atomData{'serial'}   = checkSubstr($atomStr, $strLen, 15, 5);
    $atomData{'cooX'}     = checkSubstr($atomStr, $strLen, 20, 8);
    $atomData{'cooY'}     = checkSubstr($atomStr, $strLen, 28, 8);
    $atomData{'cooZ'}     = checkSubstr($atomStr, $strLen, 36, 8);
    $atomData{'velX'}     = checkSubstr($atomStr, $strLen, 44, 8);
    $atomData{'velY'}     = checkSubstr($atomStr, $strLen, 52, 8);
    $atomData{'velZ'}     = checkSubstr($atomStr, $strLen, 60, 8);

    return \%atomData;
}



sub checkSubstr {
    my $str       = shift;
    my $strLen    = shift;
    my $start     = shift;
    my $substrLen = shift;
    my $substr    = '';

    if ($strLen >= ($start+$substrLen)) {
        $substr = substr($str, $start, $substrLen);
        $substr =~ s/\s//g;
    }
    return $substr;
}



sub writeGro {
    my $groFile    = shift;
    my $groDataRef = shift;

    $groFile .= ".gro" unless $groFile =~ /\.gro$/;

    open(GROFILE, ">$groFile") || die "ERROR: Cannot open output GRO file ($groFile): $!\n";
    writeHeader(\*GROFILE, $groDataRef);
    writeCoords(\*GROFILE, $groDataRef);
    writeFooter(\*GROFILE, $groDataRef);
    close(GROFILE);
}



sub writeHeader {
    my $fileHandle = shift;
    my $groDataRef = shift;

    print $fileHandle $$groDataRef{'title'} . "\n";
    print $fileHandle $$groDataRef{'nAtoms'} . "\n";
}



sub writeCoords {
    my $fileHandle = shift;
    my $groDataRef = shift;

    foreach (@{$$groDataRef{'atoms'}}) {
        next unless $$_{'residue'};
        printf($fileHandle "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
            ($$_{'residue'}%100000), $$_{'resName'}, $$_{'atomName'}, ($$_{'serial'}%100000), $$_{'cooX'}, $$_{'cooY'}, $$_{'cooZ'});
    }
}



sub writeFooter {
    my $fileHandle = shift;
    my $groDataRef = shift;

    printf($fileHandle "  %8.3f  %8.3f  %8.3f\n", $$groDataRef{'box'}{'cooX'}, $$groDataRef{'box'}{'cooY'}, $$groDataRef{'box'}{'cooZ'});
}

1;
