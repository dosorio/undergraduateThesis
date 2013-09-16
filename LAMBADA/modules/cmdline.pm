package Cmdline;

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
my @EXPORT = qw(addCmdlParam getCmdlParamRef parser printParamHelp);


my @cmdlParam;



sub addCmdlParam {
    push(@cmdlParam, {itype => $_[0], param => $_[1], type => $_[2], value => $_[3], defau => $_[4], descr => $_[5]});
}



sub getCmdlParamRef {
    return \@cmdlParam;
}



sub parser {
    my $argId = 0;
    my @unknownParams;

    ### Catch the command line parameters ######################################
    while ($argId<@ARGV) {
        if ($ARGV[$argId] =~ /^-(\D.*)/) {
            my $tmpParam = $1;

            my $unknown = 1; # Set the indicator for an unknown parameter as true by default.
            foreach (@cmdlParam) {
                next unless $tmpParam eq $$_{'param'} || $tmpParam eq 'no' . $$_{'param'};
                $unknown = 0; # Set the indicator for an unknown parameter back to zero.

                if ($$_{'itype'} eq 'flag') {
                    ${$$_{'value'}} = $tmpParam eq 'no' . $$_{'param'} ? 0 : 1;
                }
                elsif($$_{'itype'} eq 'array') {
                    my $tmpArgId = $argId+1;
                    while ($tmpArgId<@ARGV) {
                        if ($ARGV[$tmpArgId] =~ /^-(\D.*)/) {
                            $argId = $tmpArgId-1;
                            $tmpArgId = @ARGV;
                            next;
                        }
                        else {
                            push(@{$$_{'value'}}, $ARGV[$tmpArgId]);
                            $tmpArgId++;
                        }
                        $argId = $tmpArgId-1 if $tmpArgId == @ARGV;
                    }
                }
                else {
                    if ($ARGV[++$argId] =~ /^-(\D.*)/) {
                        ${$$_{'value'}} = $$_{'defau'};
                    }
                    else {
                        ${$$_{'value'}} = $ARGV[$argId];
                    }
                }
            }
            push(@unknownParams, $tmpParam) if $unknown;
        }
        $argId++;
    }
    ############################################################################


    ### If there are unknown parameters... #####################################
    if (@unknownParams) {
        main::printHelp();
        printParamHelp();

        print STDERR "\n-------------------------------------------------------\n";
        print STDERR "Invalid command line argument(s):\n-" . join(", -", @unknownParams) . "\n";
        print STDERR "For more information and tips for troubleshooting, please check the LAMBADA\n";
        print STDERR "website at http://code.google.com/p/lambada-align/\n";
        print STDERR "-------------------------------------------------------\n\n";
        exit;
    }
    ############################################################################

    return 1;
}



sub printParamHelp {
    my $cmdLParamRef = shift;

    printf("%5s %12s  %-11s  %s\n", "Option", "Filename", "Type", "Description");
    print ("------------------------------------------------------------\n");
    foreach (@cmdlParam) {
        next unless $$_{'type'} =~ /(In|Out)put/;
        printf("%4s  %13s  %-12s %s\n", '-' . $$_{'param'}, ${$$_{'value'}} ? ${$$_{'value'}} : $$_{'defau'}, $$_{'type'}, $$_{'descr'}) if ref($$_{'value'}) eq 'SCALAR';
        printf("%4s  %13s  %-12s %s\n", '-' . $$_{'param'}, join(" ", @{$$_{'value'}}), $$_{'type'}, $$_{'descr'}) if ref($$_{'value'}) eq 'ARRAY';
    }

    printf("\n%-12s %-6s %-6s  %s\n", "Option", "Type", "Value", "Description");
    print ("------------------------------------------------------\n");
    foreach (@cmdlParam) {
        next if $$_{'type'} =~ /(In|Out)put/;
        if ($$_{'itype'} eq 'flag') {
            printf("-%-11s %-6s %-6s  %s\n", '[no]' . $$_{'param'}, $$_{'type'}, ${$$_{'value'}} == 1 ? 'yes' : 'no',  $$_{'descr'});
        }
        else {
            printf("-%-11s %-6s %-6s  %s\n", $$_{'param'}, $$_{'type'}, ${$$_{'value'}},  $$_{'descr'});
        }
    }
    print ("\n");
}


1;