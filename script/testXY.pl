#!/usr/bin/perl 

#   ----------------------------------------------------------------------
#   written by Chen Liao, USTC
#   XY model: 
#   Hxy = -sum{j=1,N}((1+r)/2*sx(j)*sx(j+1)+(1-r)/2*sy(j)*sy(j+1)+h*sz(j))
#   Modify parameters for programme running:
#   (1) vd: virtural dimension
#   (2) ns: site number
#   (3) h:  horizontal field strength
#   (4) r:  xy anisotropy of the system
#   ----------------------------------------------------------------------

use strict;
use warnings;
use Math::Trig;

#   =======================================================================
#   @vd: array of Virtual Dims, @ns: array of site number
#   =======================================================================
my  @vd = qw\4 6 8 10 12 15 18\;
my  @ns = qw\20 40 60 80 100 120 140 160\;

#   h: horizontal field strength && paramter in XY model
my  $h = 0.96;
my  $r = 1;
#   =======================================================================
#   Plz do not modify codes below unless you know what you are doing
#   =======================================================================


























open (my $outfile, ">", "outdata_xy.txt") or die "Cound not open file outdata_xy.txt: $!";
print $outfile "NS  E1  E2  Gap\n";

#   calc gs and fes energies of XY model
my  @gs;
my  @fes;

#   periodical boundary condition
foreach my $n (@ns)
{
    my $e1 = 0.0;
    for (my $m=0; $m<$n; $m++)
    {
        my $fac1 = 2.0*pi*($m+0.5)/$n;
        $e1 -= sqrt ( ($h-cos $fac1)**2 + ($r*sin $fac1)**2 );
    }

    my $e2 = $h-1.0;
    for (my $m=1; $m<$n; $m++)
    {
        my $fac2 = 2.0*pi*$m/$n;
        $e2 -= sqrt ( ($h-cos $fac2)**2 + ($r*sin $fac2)**2 );
    }

    if ($e1 < $e2) 
    {
        push @gs, $e1;
        push @fes, $e2;
        my $gap = $e2-$e1;
        print $outfile "\n$n    $e1 $e2 $gap";     
    }
    else 
    {
        push @gs, $e2;
        push @fes, $e1;
        my $gap = $e1-$e2;
        print $outfile "\n$n    $e2 $e1 $gap";     
    }
}

print $outfile "\n\n";
#exit;

#   used to store line number
my  $ln_ns;
my  $ln_vd;
my  $ln_mth;
my  $ln_init_cfg;
my  $ln_cl;
my  $ln_j1x;
my  $ln_j1y;
my  $ln_j1z;
my  $ln_j2x;
my  $ln_j2y;
my  $ln_j2z;
my  $ln_hx;
my  $ln_hy;
my  $ln_hz;
my  $ln_tol;
my  $ln_calc_fes;
my  $ln_use_pbc; 
my  $ln_out_mag;
my  $ln_out_corr; 
my  $ln_out_cfg;

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
defined ($ln_ns) or die "Cound not find line start with NS: $!";

#   modify parameters
open (my $infile1, "<", "input.qs") or die "Cound not open file input.qs: $!";
my @lines = <$infile1>;
close $infile1;

#   modify j1x
my @new_cont_j1x = split (/\s+/, $lines[$ln_j1x-1]);
my $val1 = (1.0+$r)/2.0;
$new_cont_j1x[1] = "$val1\n";
$lines[$ln_j1x-1] = join ("      ", @new_cont_j1x);

#   modify j1y
my @new_cont_j1y = split (/\s+/, $lines[$ln_j1y-1]);
my $val2 = (1.0-$r)/2.0;
$new_cont_j1y[1] = "$val2\n";
$lines[$ln_j1y-1] = join ("      ", @new_cont_j1y);

#   modify hz
my @new_cont_hz = split (/\s+/, $lines[$ln_hz-1]);
$new_cont_hz[1] = "$h\n";
$lines[$ln_hz-1] = join ("      ", @new_cont_hz);

#   modify use_pbc 
my @new_cont_pbc = split (/\s+/, $lines[$ln_use_pbc-1]);
$new_cont_pbc[1] = "1\n";
$lines[$ln_use_pbc-1] = join ("      ", @new_cont_pbc);

open (my $outfile1, ">", "input.qs") or die "Cound not open file input.qs: $!";
print {$outfile1} @lines;
close $outfile1;

#   the best way to modify contents of a file is to 
#   fetch each line into an array and modfify its contents
#   and then write to the file

foreach my $curr_v (@vd)
{
    foreach my $curr_n (@ns)
    {
        open (my $infile2, "<", "input.qs") or die "Cound not open file input.qs: $!";
        my @lines = <$infile2>;
        close $infile2;

        #  modify VD
        my @new_cont_xd = split (/\s+/, $lines[$ln_vd-1]);
        my $num = scalar(@new_cont_xd);
        $new_cont_xd[1] = "$curr_v\n";
        $lines[$ln_vd-1] = join ("      ", @new_cont_xd);

        #   modify NS
        my @new_cont_ns = split (/\s+/, $lines[$ln_ns-1]);
        $new_cont_ns[1] = "$curr_n\n";
        $lines[$ln_ns-1] = join ("      ", @new_cont_ns);

        #   output to file
        open (my $outfile2, ">", "input.qs") or die "Cound not open file input.qs: $!";
        print {$outfile2} @lines;
        close $outfile2;

        system ("./qsAnnni");

        #   fetch data
        system ("cat output_res.qs >> outdata_xy.txt");
    }
}

close $outfile;
