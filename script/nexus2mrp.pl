#!/usr/bin/perl
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use Bio::Phylo::IO 'parse';

# process command line arguments
my ( $infile, $outfile, $wdir );
GetOptions(
	'infile=s'  => \$infile,
	'outfile=s' => \$outfile,
	'dir=s'     => \$wdir,		
);

# parse nexus
my $forest = parse(
	'-format' => 'nexus',
	'-file'   => $infile,
)->[-1];

# make base name
my ( $vol, $dir, $file ) = File::Spec->splitpath( $infile );
$file =~ s/\..+$//;

# make outfile name
if ( not $outfile ) {
	my $odir = $wdir || $dir;
	$outfile = File::Spec->catfile( $odir, $file . '.dat' );
}

# write result
open my $fh, '>', $outfile or die $!;
$forest->make_matrix->visit(sub{
	my $row  = shift;
	my $char = $row->get_char;
	my $name = $row->get_name;
	print $fh $file, "\t", $name, "\t", $char, "\n";
});