#!/usr/bin/perl

=head1 Name

 phylo_tree.pl - Using TreeBeST or PhyML or MrBayes to 
 reconstruct phylogenetic tree

=head1 Usage
 
 perl phylo_tree.pl <infile> [options]
  <infile>         input phylip format(sequential) or mfa format sequence file;
  -format <str>    file formats: phylip(default); mfa;
  -type <str>      method: nj(default); ml; bayes;
  -outdir <str>    save results to this directory, default ./;
  -clean           delete temporary files, default not;
  -help            show help information;
 
 TreeBeST's relevant parameters(nj):
  -d <str>         data type: nt, nucleotide(default); aa, amino acid;
  -t <str>         model for codon & nt: ntmm; dn; ds; dm(default); 
                   model for aa: mm(default); jtt; kimura;
  -b <num>         bootstrap times, default 100; 

 PhyML's relevant parameters(ml):
  -d <str>         data type: nt, nucleotide(default); aa, amino acid;
  -m <str>         models for nt: HKY85(default); JC69; K80; F81; F84; TN93; GTR;
                   models for aa: WAG(default); JTT; MtREV; Dayhoff;
  -b <num>         set bootstrap(int>0) or aLRT method(-1,-2,-4), default -4;
  -rates <str>     rates across sites: gamma(default); invgamma; equal; propinv;

 MrBayes's relevant parameters(bayes):
  -data <str>      data type: dna(default); protein;
  -nst <num>       nucleotide models: 6, GTR(default); 2, HKY; 1, JC;
  -aamodel <str>   amino-acid models: Poisson; Jones(default); Dayhoff; WAG; BLOSUM;
  -rates <str>     rates across sites: gamma(default); invgamma; equal; propinv;
  -ngen <num>      generation number, default 100000;
  -samplefreq <num>    sample one time for every 100 (default) generations;
  -mpirun <num>    run MrBayes with mpirun using multi CPUs or not, default 1;

=head1 Example 

 perl phylo_tree.pl <in.phy> -type nj -d nt
 perl phylo_tree.pl <in.phy> -type ml -d nt
 perl phylo_tree.pl <in.phy> -type bayes -data dna

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin";
use lib "$Bin/../lib";
use Cwd;

my ( $type, $outdir, $format, $mpirun ) = ( "nj", "./", "phylip", 1 );
my ( $help, $clean );
my ( $aamodel, $data, $rates, $nst, $ngen, $samplefreq);
my ( $d, $m, $b, $t );
GetOptions(
    'format:s' => \$format,
    'type:s'   => \$type,
    'outdir:s' => \$outdir,
    'clean'    => \$clean,
    'help'     => \$help,

    'data:s'       => \$data,
    'aamodel:s'    => \$aamodel,
    'nst:i'        => \$nst,
    'rates:s'      => \$rates,
    'ngen:i'       => \$ngen,
    'samplefreq:i' => \$samplefreq,
    'mpirun:i'     => \$mpirun,

    'd:s' => \$d,
    't:s' => \$t,
    'b:i' => \$b,

    'm:s' => \$m,
);
die `pod2text $0` if ( @ARGV < 1 || $help );
my $infile = shift;
mkdir $outdir if ( !-e $outdir );
$format ||= "phylip";

if ( $type =~ /bayes/ ) {
    $data       ||= "dna";
    $rates      ||= "gamma";
    $nst        ||= 6;
    $ngen       ||= 100000;
    $samplefreq ||= 100;
    $aamodel    ||= "Jones";
    $mpirun;
    my $burnin = ( $ngen / $samplefreq ) * 0.25;

    if ( $format =~ /mfa/ ) {
        my $number = `grep ">" $infile -c`;
        chomp $number;
        my ( %hash, $name, $sequence, $length );
        open( IN, "$infile" ) or die "Can't open: $infile\n";
        open( OUT, ">$outdir/temp.phy" )
          or die "Can't open: $outdir/temp.phy\n";
        while (<IN>) {
            chomp;
            if ( $_ =~ s/^\>(\w+)/$1/ ) {
                $name = $1;
            }
            else {
                $sequence = $_;
                $hash{$name} .= $sequence;
            }
            $length = length( $hash{$name} );
        }
        print OUT "    $number         $length\n";
        foreach $name ( keys %hash ) {
            print OUT "$name        $hash{$name}\n";
        }
        close IN;
        close OUT;
    }
    elsif ( $format =~ /phylip/ ) {
        `cp $infile $outdir/temp.phy`;
    }
    open( IN,  "$outdir/temp.phy" )   or die "Can't open: $outdir/temp.phy\n";
    open( OUT, ">$outdir/input.nex" ) or die "Can't open: $outdir/input.nex\n";
    my $firstline = <IN>;
    my ( $species_number, $sequence_length );
    if ( $firstline =~ /^\s*(\d+)\s+(\d+)/ ) {
        $species_number = $1, $sequence_length = $2;
        print OUT
"#NEXUS\n\nbegin data;\n\tdimensions ntax =$species_number nchar=$sequence_length;";
    }
    if ( $sequence_length < 15000 ) {
        print OUT
"\n\tformat datatype=$data interleave =no missing=? gap=-;\n\tmatrix\n";
        while ( my $eachline = <IN> ) {
            if ( $eachline =~ /^\w+\s+[\w-]/ ) {
                print OUT "\t$eachline";
            }
        }
    }
    else {
        my ( @phylip, @name, @sequence, @aim, $change );
        @phylip = <IN>;
        foreach (@phylip) {
            my @temp = split /\s+/, $_;
            push @name,     $temp[0];
            push @sequence, $temp[1];
        }
        for ( my $i = 0 ; $i < $sequence_length ; $i += 15000 ) {
            for ( my $k = 0 ; $k < scalar(@sequence) ; $k++ ) {
                chomp $sequence[$k];
                my $change = substr( $sequence[$k], $i, 15000 );
                push @aim, "\t$name[$k]\t$change\n";
            }
        }
        print OUT
"\n\tformat datatype=$data interleave =yes missing=? gap=-;\n\tmatrix\n";
        foreach (@aim) {
            chomp;
            print OUT "$_\n";
        }
        close IN;
    }
    print OUT "\t;\nend;\n";
    print OUT "begin mrbayes;\n";
    print OUT "\tset autoclose=yes nowarn=yes;\n";
    if ( lc($data) eq "dna" ) {
        print OUT "\tlset nst=$nst rates=$rates;\n";
    }
    elsif ( lc($data) eq "protein" ) {
        print OUT "\tprset aamodelpr=fixed($aamodel);\n\tlset rates=$rates;\n";
    }
    print OUT
"\tmcmc ngen=$ngen printfreq=$samplefreq samplefreq=$samplefreq nchains=4 savebrlens=yes;\n";
    print OUT
      "\tsumt filename=$outdir/input.nex contype=allcompat burnin=$burnin;\n";
    print OUT "\tsump filename=$outdir/input.nex burnin=$burnin printtofile=yes;\n";
    print OUT "\tquit;\nend;\n";
    close OUT;

    if ($mpirun > 1) {
        system "mpirun -c $mpirun $Bin/software/mb_mpi  -i $outdir/input.nex ";

#system "mpirun -np 2 -mca mpi_show_handle_leaks 1 $Bin/software/mb  -i $outdir/input.nex ";
    }
    else {
        system "$Bin/software/mb  -i $outdir/input.nex ";
    }

    print "Note: generation number must be large enough until \n",
      "'Average standard deviation of split frequencies' less than 0.01\n";

    open( IN, "$outdir/input.nex.con" )
      or die "Can't open: $outdir/input.nex.con\n";
    my $treestr;
    while (<IN>) {
        chomp;
        if ( $_ =~ /^\s.+=/ ) {
            my @arr = split /=/, $_;
            $treestr = $arr[1];
            last;
        }
    }
    close IN;
    open( OUT, ">$outdir/tree.newick" )
      or die "Can't open: $outdir/tree.newick\n";
    my $str = $treestr;
    $str =~ s/\)\d.\d+\:/\)\:/g;
    print OUT "$str\n";
    close OUT;
    ## To get B value ##
    open( OUT, ">$outdir/tree.nhx" ) or die "Can't open: $outdir/tree.nhx\n";
    my $str = $treestr;
    $str =~ s/\)([\d|.]+)\:([\d|.]+)/res($1, $2)/ge;
    print OUT $str;
    close OUT;

    sub res {
        my ( $a, $b ) = @_;
        my $c = $a * 100;
        return "\)\:$b\[&&NHX:B=$c\]";
    }
    `rm $outdir/temp.phy`;
    `rm  $outdir/tree.newick $outdir/input.nex*` if ($clean);

}
elsif ( $type =~ /nj/ ) {

    $d ||= "nt";
    if    ( $d eq "nt" ) { $t ||= "dm"; }
    elsif ( $d eq "aa" ) { $t ||= "mm"; }
    $b ||= "100";
    if ( $format =~ /phylip/ ) {
        open( IN, "$infile" ) or die "Can't open: $infile\n";
        <IN>;
        open( OUT, ">$outdir/in.mfa" ) or die "Can't open: $outdir/in.mfa\n";
        while (<IN>) {
            chomp;
            my @arr1 = split /\s+/;
            print OUT ">$arr1[0]\n$arr1[1]\n";
        }
        close IN;
        close OUT;
    }
    elsif ( $format =~ /mfa/ ) {
        `cp  $infile $outdir/in.mfa`;
    }
    print "datatype = $d\nmodel = $t\nbootstrap = $b\n";
    ###call njtree##
    system
"$Bin/software/treebest nj -t $t $outdir/in.mfa -b $b 1> $outdir/tree.nhx 2> $outdir/err.log";
    open( IN, "$outdir/tree.nhx" ) or die "Can't open: $outdir/tree.nhx\n";
    my @arr2 = <IN>;
    close IN;
    my $treestr = join( "", @arr2 );
    $treestr =~ s/[\n\r]+//g;

    open( OUT, ">$outdir/tree.newick" )
      or die "Can't open: $outdir/tree.newick\n";
    my $str = $treestr;
    $str =~ s/\[.+?\]//g;
    print OUT "$str\n";
    close OUT;

    open( OUT, ">$outdir/tree.nhx" ) or die "Can't open: $outdir/tree.nhx\n";
    print OUT "$treestr\n";
    close OUT;
    `rm $outdir/tree.newick $outdir/in.mfa` if ($clean);

}
elsif ( $type =~ /ml/ ) {
    $d     ||= "nt";
    $b     ||= -4;
    $rates ||= "gamma";
    if ( $format =~ /mfa/ ) {
        my $number = `grep ">" $infile -c`;
        chomp $number;
        my ( %hash, $name, $sequence, $length );
        open( IN, "$infile" ) or die "Can't open: $infile\n";
        open( OUT, ">$outdir/temp.phy" )
          or die "Can't open: $outdir/temp.phy\n";
        while (<IN>) {
            chomp;
            if ( $_ =~ s/\>(\w+)/$1/ ) {
                $name = $1;
            }
            else {
                $sequence = $_;
                $hash{$name} .= $sequence;
            }
            $length = length( $hash{$name} );
        }
        print OUT "    $number     $length\n";
        foreach $name ( keys %hash ) {
            print OUT "$name        $hash{$name}\n";
        }
        close IN;
        close OUT;
    }
    elsif ( $format =~ /phylip/ ) {
        `cp $infile $outdir/temp.phy`;
    }

    my $cwd = Cwd::getcwd();
    chdir "$outdir/";

    if ( $d eq "nt" ) {
        $m ||= "HKY85";
        my $same = "$Bin/software/phyml -i temp.phy -d $d -b $b -m $m ";
        if    ( $rates =~ /equal/ )    { system "$same -t e"; }
        elsif ( $rates =~ /invgamma/ ) { system "$same -a e -c 4 -v e -t e"; }
        elsif ( $rates =~ /gamma/ )    { system "$same -a e -c 4 -t e"; }
        elsif ( $rates =~ /propinv/ )  { system "$same -v e -t e"; }
    }
    elsif ( $d eq "aa" ) {
        $m ||= "WAG";
        my $same = "$Bin/software/phyml -i temp.phy -d $d -b $b -m $m ";
        if    ( $rates =~ /equal/ )    { system "$same"; }
        elsif ( $rates =~ /invgamma/ ) { system "$same -a e -c 4 -v e"; }
        elsif ( $rates =~ /gamma/ )    { system "$same -a e -c 4"; }
        elsif ( $rates =~ /propinv/ )  { system "$same -v e"; }

    }
    chdir "$cwd";
    open( IN, "$outdir/temp.phy_phyml_tree.txt" )
      or die "Can't open: $outdir/temp.phy_phyml_tree.txt\n";
    open( OUT, ">$outdir/tree.newick" )
      or die "Can't open: $outdir/tree.newick\n";
    my $overall = <IN>;
    my $treestr = $overall;
    if ( $b < 0 ) {
        $treestr =~ s/\)\d.\d+\:/\)\:/g;
        print OUT "$treestr";
    }
    else {
        $treestr =~ s/\)\d+\:/\)\:/g;
        print OUT "$treestr";
    }
    close OUT;

    ##To get B value##
    open( OUT, ">$outdir/tree.nhx" ) or die "Can't open: $outdir/tree.nhx\n";
    my $treestr1 = $overall;
    chomp $treestr1;
    if ( $b > 0 ) {
        $treestr1 =~ s/\)(\d+)\:([\d|.]+)/\)\:$2\[&&NHX:B=$1\]/g;
        print OUT"$treestr1\n";
    }
    else {
        $treestr1 =~ s/\)([\d|.]+)\:([\d|.]+)/res($1, $2)/ge;
        print OUT "$treestr1";

        sub res {
            my ( $a, $b ) = @_;
            my $c = $a * 100;
            return "\)\:$b\[&&NHX:B=$c\]";
        }
    }
    close IN;
    close OUT;
    `rm $outdir/tree.newick $outdir/temp.phy $outdir/temp.phy_phyml_*`
      if ($clean);
}
system "perl $Bin/draw_tree.pl $outdir/tree.nhx > $outdir/tree.svg";
