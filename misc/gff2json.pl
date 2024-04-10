#!/usr/bin/env perl

use warnings;
use strict;
use lib "/data/molgentools/lib";
use anne_files ;
use anne_misc ;
use anne_genomics ;
use JSON;


# ------------------------------------------------------------ parameters ----------------------------------------------------------------------


my $query ;
my $type = 'gene|CDS|RNA' ;

my $usage = "/data/ppp/gff2json.pl
	-gff GFF3 file
	-type type to include [default= $type]
	Result will be a json file

";



&parseparam() ;


# ------------------------------------------------------------ main ----------------------------------------------------------------------

	
	my %gff = anne_genomics::gff_2_hash_v3($query, $type) ;
	my $outputfile = $query ;
	$outputfile =~ s/gff$/json/ ;
	anne_files::write_lines($outputfile,(encode_json(\%gff))) ;
	
	my @genearray ;
	foreach my $genome (keys %gff) {
		foreach my $type (keys $gff{$genome}) {
			foreach my $ID (keys $gff{$genome}{$type}) {
				my $gene = "\n{";
				$gene .= '"Name": "'.$gff{$genome}{$type}{$ID}{Name}."\",\n";
				$gene .= '"start": '.$gff{$genome}{$type}{$ID}{start}.",\n";
				$gene .= '"end": '.$gff{$genome}{$type}{$ID}{end}.",\n";
				$gene .= '"strand": "'.$gff{$genome}{$type}{$ID}{strand}."\"\n";
				$gene .= "}";
				push @genearray, $gene ;
			}
		}
	}
	my $gene_json = join ',',@genearray ;	
	$outputfile =~ s/json$/2\.json/ ;
	anne_files::write_lines($outputfile,"[ $gene_json ]") ;

# ------------------------------------------------------------ functions ----------------------------------------------------------------------


sub parseparam {
    my $var ;
    my @arg = @ARGV ;

    while(@arg) {
        $var = shift(@arg) ;
		die $usage if ($var eq '-h' or $var eq '--help') ;
		$query 	= shift(@arg) if($var eq '-gff') ;
		$type 	= shift(@arg) if($var eq '-type') ;
    }
    die $usage if (!$query) ;
}
