=head1 LICENSE

Copyright [2026] Romualdo Morandi Filho

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

 romualdofarm@gmail.com
    
=cut

=head1 NAME

 ABraOM

=head1 SYNOPSIS

 mv ABraOM.pm ~/.vep/Plugins

 Annotate only the default column abraom_freq (from the original column Frequencies):
 ./vep -i variations.vcf --plugin ABraOM,file=/path/to/abraom/SABE1171.Abraom_final.tsv.gz

 Annotate any specified columns (e.g. Frequencies, Gene.refGene, PredictedFunc.refGene) that will be added with the prefix "abraom_":
 ./vep -i variations.vcf --plugin ABraOM,file=/path/to/abraom/SABE1171.Abraom_final.tsv.gz,cols=Frequencies,Gene.refGene,PredictedFunc.refGene

=head1 DESCRIPTION

 An Ensembl VEP plugin that adds transcript-level population frequency annotations from the ABraOM database.

 The ABraOM repository contains genomic variants obtained with whole-exome and whole-genome sequencing from SABE,
 a census-based sample of elderly individuals from São Paulo, Brazil's largest city.
 
 Please cite the ABraOM publication alongside the VEP if you use this resource:
 http://onlinelibrary.wiley.com/doi/10.1002/humu.23220/full

 The database file can be downloaded from
 https://abraom.ib.usp.br/download/

 This file can be tabix-processed by:
 tar -xvzf SABE1171.Abraom.clean.tar.gz
 cat SABE1171.Abraom.clean.tar.gz | (head -n 1 && tail -n +2  | sort -t$'\t' -k 1,1 -k 2,2n ) > SABE1171.Abraom_sorted.tsv
 sed '1s/.*/#&/' SABE1171.Abraom_sorted.tsv > SABE1171.Abraom_final.tsv
 bgzip SABE1171.Abraom_final.tsv
 tabix -f -s 1 -b 2 -e 2 SABE1171.Abraom_final.tsv.gz

=cut

package abraom;

use strict;
use warnings;
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);
use parent qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub new {
  my ($class, @params) = @_;
  my $self = $class->SUPER::new(@params);
  bless $self, $class;

  # parameters
  my ($file, $cols) = (undef, undef);
  foreach my $p (@params) {
    if ($p =~ /^file=(.+)$/) {
      $file = $1;
    } elsif ($p =~ /^cols=(.+)$/) {
      $cols = $1;
    }
  }
  die("You need an ABraOM .gz with 'file='") unless $file;

  my @columns = $cols ? split(/,/, $cols) : ('Frequencies');
  $self->{columns} = \@columns;

  # inicializa tabix via BaseVepTabixPlugin
  $self->add_file($file);

  return $self;
}

sub get_header_info {
  my ($self) = @_;
  my %header;
  foreach my $c (@{$self->{columns}}) {
    my $tag = $c eq 'Frequencies' ? 'abraom_freq' : "abraom_$c";
    $header{$tag} = "Column $c from ABraOM database";
  }
  return \%header;
}

sub parse_data {
  my ($self, $line, $file) = @_;
  my @f = split(/\t/, $line);
  
  # Mapeamento de índices das colunas do arquivo ABraOM
  my %col_index = (
    'Chr' => 0,
    'Start' => 1,
    'Ref' => 2,
    'Alt' => 3,
    'PredictedFunc.refGene' => 4,
    'Gene.refGene' => 5,
    'PredConsequence.refGene' => 6,
    'avsnp150' => 7,
    'FILTER' => 8,
    'CEGH_Filter' => 9,
    'HomozygousALT_count' => 10,
    'Hemizygous_count' => 11,
    'Allele_number' => 12,
    'Allele_ALT_count' => 13,
    'Frequencies' => 14,
    'Cohort' => 15,
  );
  
  # Extract data for the requested columns
  my %result = ();
  foreach my $c (@{$self->{columns}}) {
    my $tag = $c eq 'Frequencies' ? 'abraom_freq' : "abraom_$c";
    my $i = $col_index{$c};
    $result{$tag} = $f[$i] if defined $i;
  }
  
  return {
    ref => $f[2],
    alt => $f[3],
    pos => $f[1],
    result => \%result
  };
}

sub get_start {
  my ($self, $data) = @_;
  return $data->{pos};
}

sub get_end {
  my ($self, $data) = @_;
  return $data->{pos};
}

sub run {
  my ($self, $tva) = @_;
  my $vf = $tva->variation_feature;
  my $chr = $vf->seq_region_name;
  my $pos = $vf->start;
  my $ref = $vf->ref_allele_string;
  my $alt = $tva->allele_string;
  
  # Extract only the alternative allele if it is in "REF/ALT" format
  $alt =~ s/^[^\/]*\/// if $alt =~ /\//;

  # get_data() returns an arrayref of parsed data via parse_data()
  my @data = @{$self->get_data($chr, $pos, $pos)};
  return {} unless @data;

  foreach my $entry (@data) {
    # Comparison of alleles between VEP variant and ABraOM variant
    my $var_a = {
      ref    => $ref,
      alts   => [$alt],
      pos    => $pos,
      strand => $vf->strand,
    };
    my $var_b = {
      ref  => $entry->{ref},
      alts => [$entry->{alt}],
      pos  => $entry->{pos},
    };
    my @matched = @{get_matched_variant_alleles($var_a, $var_b)};
    next unless @matched;
    
    # Return the results
    return $entry->{result};
  }

  return {};
}
1
