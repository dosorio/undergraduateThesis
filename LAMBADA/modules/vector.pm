package Vector;

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
my @EXPORT = qw(vSub vLen vNorm vXprod vDotprod);


sub vSub {
    my %vecDiff = ('cooX' => $_[0]{'cooX'} - $_[1]{'cooX'},
                   'cooY' => $_[0]{'cooY'} - $_[1]{'cooY'},
                   'cooZ' => $_[0]{'cooZ'} - $_[1]{'cooZ'});
    return %vecDiff;
}



sub vLen {
    return sqrt($_[0]{'cooX'}*$_[0]{'cooX'} + $_[0]{'cooY'}*$_[0]{'cooY'} + $_[0]{'cooZ'}*$_[0]{'cooZ'});
}



sub vNorm {
    my $vecLen  = vLen($_[0]);
    my %vecNorm = ('cooX' => $_[0]{'cooX'}/$vecLen,
                   'cooY' => $_[0]{'cooY'}/$vecLen,
                   'cooZ' => $_[0]{'cooZ'}/$vecLen);
    return %vecNorm;
}



sub vXprod {
    my %vecXprod = ('cooX' => $_[0]{'cooY'} * $_[1]{'cooZ'} - $_[0]{'cooZ'} * $_[1]{'cooY'},
                    'cooY' => $_[0]{'cooZ'} * $_[1]{'cooX'} - $_[0]{'cooX'} * $_[1]{'cooZ'},
                    'cooZ' => $_[0]{'cooX'} * $_[1]{'cooY'} - $_[0]{'cooY'} * $_[1]{'cooX'});
    return %vecXprod;
}



sub vDotprod {
    return ($_[0]{'cooX'} * $_[1]{'cooX'} + $_[0]{'cooY'} * $_[1]{'cooY'} + $_[0]{'cooZ'} * $_[1]{'cooZ'});
}

1;