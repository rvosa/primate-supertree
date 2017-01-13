#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Bio::Align::DNAStatistics;
use Bio::Tree::DistanceFactory;
use Bio::Phylo::Factory;
use Bio::Phylo::Forest::Tree;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::IO qw'parse_tree unparse';

# process command line arguments
my $infile;
my $verbosity = WARN;
GetOptions(
	'infile=s' => \$infile,
	'verbose+' => \$verbosity,
);

# instantiate helper objects
my $fac = Bio::Phylo::Factory->new;
my $dna = Bio::Align::DNAStatistics->new;
my $njf = Bio::Tree::DistanceFactory->new( '-method' => "NJ" );
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);

# start reading the file
my $block;
my %data;
open my $fh, '<', $infile or die $!;
while(<$fh>) {
	chomp;
	my ( $tb, $name, $char ) = split /\t/, $_;
	if ( not $block ) {
		$block = $tb;
		$data{$name} = $char;
	}
	else {
		if ( $block ne $tb ) {
			$log->info("going to analyze block $block");
			analyze_block( %data );
			%data = ( $name => $char );
			$block = $tb;
		}
		else {
			$data{$name} = $char;
		}	
	}
}
$log->info("going to analyze last block $block");
$log->debug(Dumper(\%data));
analyze_block( %data );

sub analyze_block {
	my %rows = @_;
	
	# make outgroup if needed
	my ($row) = values %rows;
	my $nchar = length $row;
	my ($out) = grep { /outgroup/ } keys %rows;
	if ( not $out ) {
		$out = 'myoutgroup';
		$rows{$out} = '0' x $nchar;
		$log->info("creating outgroup with $nchar characters");
	}
	
	# find the partitions
	my $previous = 0;
	for my $i ( 0 .. ( $nchar - 1 ) ) {
		my %c;
		ROW: for my $row ( keys %rows ) {
			next ROW if $row eq $out;
			my $c = substr $rows{$row}, $i, 1;
			$c{$c} = 1;
		}
		
		# first column is also all 1's but we don't care
		if ( $i != 0 ) {
		
			# found next column with all 1's
			if ( 1 == scalar keys %c ) {
				$log->info("going to analyze matrix $previous..$i");
				analyze_matrix( $previous, $i, %rows );			
				$previous = $i;		
			}
			
			# reached end of data
			elsif ( $i == ( $nchar - 1 ) ) {
				$log->info("going to analyze final matrix $previous..$nchar");
				analyze_matrix( $previous, $nchar, %rows );
			}
		}
	}
}

sub analyze_matrix {
	my ( $start, $end, %data ) = @_;
	
	# populate matrix
	my $m = $fac->create_matrix( '-type' => 'dna' );
	for my $row ( keys %data ) {
		my $char = substr $data{$row}, $start, ( $end - $start );
		$char =~ tr/012/AC?/;
		$m->insert(
			$fac->create_datum(
				'-type' => 'dna',
				'-name' => $row,
				'-char' => $char,
			)
		);
		$log->debug("added $row => $char");
	}
	
	# do NJ tree
	my $dist = $dna->distance(
		'-align'  => $m,
		'-method' => 'Uncorrected',
	);
	my $tree = Bio::Phylo::Forest::Tree->new_from_bioperl($njf->make_tree($dist));
	$log->info($tree->to_newick);
	make_tree_output($tree);
}

sub make_tree_output {
	my $tree = shift;
	$log->info("going to make distance matrix");
	my ($out) = grep { $_->get_name =~ /outgroup/ } @{ $tree->get_terminals };
	$out->set_root_below;
	$tree->prune_tips([$out]);
	$tree->visit(sub{
		my $n = shift;
		if ( $n->is_root ) {
			$n->set_branch_length();
		}
		else {
			$n->set_branch_length(1);
		}
	});
	print $tree->to_newick, "\n";
}

