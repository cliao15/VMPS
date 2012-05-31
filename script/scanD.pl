#!/usr/bin/perl 
#   written by Chen Liao, USTC

use strict;
use warnings;

#   =======================================================================
#   @v: array of Virtual Dims, @h: array of horizontal field strength
#   =======================================================================

#my @v = qw\2 3 5 8 10 12 15\;
my @v = qw\2 3\;    #   VD
my @h = qw\0.95 0.96 0.97 0.98 0.99 1.00 1.01 1.02 1.03 1.04 1.05\; #   h

#   =======================================================================
#   Plz do not modify codes below unless you know what you are doing
#   =======================================================================

#   used to store line number
my $ln_ns;
my $ln_vd;
my $ln_mth;
my $ln_init_cfg;
my $ln_cl;
my $ln_j1x;
my $ln_j1y;
my $ln_j1z;
my $ln_j2x;
my $ln_j2y;
my $ln_j2z;
my $ln_hx;
my $ln_hy;
my $ln_hz;
my $ln_tol;
my $ln_calc_fes;
my $ln_use_pbc; 
my $ln_out_mag;
my $ln_out_corr; 
my $ln_out_cfg;

#   find line number of 
my $inputf;
open ($inputf, "< input.qs") or die "Cound not open file input.qs: $!";
while (my $line = <$inputf>)
{
    if ($line =~ /^NS/) {$ln_ns = $.;}
    if ($line =~ /^VD/) {$ln_vd = $.;}
    if ($line =~ /^MTH/) {$ln_mth = $.;}
    if ($line =~ /^INIT_CFG/) {$ln_init_cfg = $.;}
    if ($line =~ /^CL/) {$ln_cl = $.;}
    if ($line =~ /^J1X/) {$ln_j1x = $.;}
    if ($line =~ /^J1Y/) {$ln_j1y = $.;}
    if ($line =~ /^J1Z/) {$ln_j1z = $.;}
    if ($line =~ /^J2X/) {$ln_j2x = $.;}
    if ($line =~ /^J2Y/) {$ln_j2y = $.;}
    if ($line =~ /^J2Z/) {$ln_j2z = $.;}
    if ($line =~ /^HX/) {$ln_hx = $.;}
    if ($line =~ /^HY/) {$ln_hy = $.;}
    if ($line =~ /^HZ/) {$ln_hz = $.;}
    if ($line =~ /^TOL/) {$ln_tol = $.;}
    if ($line =~ /^CALC_FES/) {$ln_calc_fes = $.;}
    if ($line =~ /^USE_PBC/) {$ln_use_pbc = $.;}
    if ($line =~ /^OUT_MAG/) {$ln_out_mag = $.;}
    if ($line =~ /^OUT_CORR/) {$ln_out_corr = $.;}
    if ($line =~ /^OUT_CFG/) {$ln_out_cfg = $.;}
}

close ($inputf);

#   defined VD and HZ
defined ($ln_vd) or die "Cound not find line start with VD: $!";
defined ($ln_hz) or die "Cound not find line start with HZ: $!";

#   rm outdata.file
if (-e "outdata.txt") {system ("rm outdata.txt");}

#   the best way to modify contents of a file is to 
#   fetch each line into an array and modfify its contents
#   and then write to the file

foreach my $curr_v (@v)
{
    foreach my $curr_h (@h)
    {
        open (my $infile, "<", "input.qs") or die "Cound not open file input.qs: $!";
        my @lines = <$infile>;
        close $infile;

        #  modify VD
        my @new_cont_xd = split (/\s+/, $lines[$ln_vd-1]);
        my $num = scalar(@new_cont_xd);
        $new_cont_xd[1] = "$curr_v\n";
        $lines[$ln_vd-1] = join ("      ", @new_cont_xd);
#        print "$lines[$ln_vd-1]\n";

        #   modify HZ
        my @new_cont_hz = split (/\s+/, $lines[$ln_hz-1]);
        $new_cont_hz[1] = "$curr_h\n";
        $lines[$ln_hz-1] = join ("      ", @new_cont_hz);
#        print "$lines[$ln_hz-1]\n";

        #   output to file
        open (my $outfile, ">", "input.qs") or die "Cound not open file input.qs: $!";
        print {$outfile} @lines;

        system ("./qsAnnni");

        #   fetch data
        system ("cat output_res.qs >> outdata.txt");

        close $outfile;
    }
}
