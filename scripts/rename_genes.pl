#!/usr/bin/env perl
#
# rename_genes.pl - assign systematic gene/mRNA IDs to a PASA/EVM GFF3
#
# Naming scheme (default):
#
#     gene   Cla05AG00160
#            ^^^ ^^ ^ ^^^^^
#            |   |  | |
#            |   |  | +-- 5-digit serial, numbered along the chromosome in
#            |   |  |     steps of 10 so new genes can be inserted later
#            |   |  +---- literal 'G'
#            |   +------- haplotype A / B / C
#            +----------- 2-digit chromosome, zero padded
#     ^^^ species prefix (Cla = Clausena lansium)
#
#     mRNA   Cla05AG00160.1   .2 .3 ...   isoforms, longest CDS first
#     exon   Cla05AG00160.1.exon1
#     CDS    Cla05AG00160.1.cds1
#     UTR    Cla05AG00160.1.utr5p1 / .utr3p1
#
# The old identifier is kept as  prev_id=  on every gene and mRNA, so the new
# annotation can always be traced back to the EVM/PASA models.
#
# Usage:
#   perl rename_genes.pl -i CL_pasa.gff3 -o CL.renamed.gff3 -p Cla
#   perl rename_genes.pl -i in.gff3 -o out.gff3 -p Cla --step 10 --digits 5 \
#                        -m id_map.tsv
#
# Options:
#   -i, --in      FILE   input GFF3                         [required]
#   -o, --out     FILE   output GFF3                        [required]
#   -m, --map     FILE   old->new mapping table             [<out>.idmap.tsv]
#   -p, --prefix  STR    species prefix                     [Cla]
#       --step    INT    increment between genes            [10]
#       --digits  INT    digits in the serial               [5]
#       --source  STR    force column 2 to this value       [keep as is]
#       --keep-name      keep the original Name= attribute  [default: drop it]
#   -h, --help
#
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my ($in, $out, $mapfile, $source, $keep_name, $help);
my $prefix = 'Cla';
my $step   = 10;
my $digits = 5;

GetOptions(
    'i|in=s'     => \$in,
    'o|out=s'    => \$out,
    'm|map=s'    => \$mapfile,
    'p|prefix=s' => \$prefix,
    'step=i'     => \$step,
    'digits=i'   => \$digits,
    'source=s'   => \$source,
    'keep-name'  => \$keep_name,
    'h|help'     => \$help,
) or die "Bad options; try --help\n";
usage() if $help;
usage('-i and -o are required.') unless $in && $out;
$mapfile = "$out.idmap.tsv" unless defined $mapfile;

##------------------------------------------------------------------
## 1. read the GFF3
##    Comment lines (# PASA_UPDATE, #PROT, ...) are dropped: #PROT lines
##    carry the old protein sequence and would be stale after renaming.
##------------------------------------------------------------------
my (%gene, %mrna, %child, @gene_order);

my $IN;
if ($in =~ /\.gz$/) { open $IN, '-|', "gzip -dc \Q$in\E" or die "$in: $!\n" }
else                { open $IN, '<', $in                 or die "$in: $!\n" }

my $nline = 0;
while (my $l = <$IN>) {
    $l =~ s/\r?\n\z//;
    next if $l =~ /^#/ || $l =~ /^\s*$/;
    my @f = split /\t/, $l;
    next if @f < 9;
    $nline++;

    my ($id)     = $f[8] =~ /\bID=([^;]+)/;
    my ($parent) = $f[8] =~ /\bParent=([^;]+)/;

    if ($f[2] eq 'gene') {
        die "Line $.: gene without ID\n" unless defined $id;
        $gene{$id} = { f => [@f], mrnas => [] };
        push @gene_order, $id;
    }
    elsif ($f[2] eq 'mRNA' || $f[2] eq 'transcript') {
        die "Line $.: mRNA without ID\n" unless defined $id;
        $mrna{$id} = { f => [@f], parent => $parent, feats => [] };
        push @{ $gene{$parent}{mrnas} }, $id if defined $parent && $gene{$parent};
    }
    else {
        next unless defined $parent;
        push @{ $child{$parent} }, [@f];
    }
}
close $IN;
die "No features found in $in\n" unless @gene_order;
warn "read $nline feature lines: ", scalar(@gene_order), " genes, ",
     scalar(keys %mrna), " mRNAs\n";

##------------------------------------------------------------------
## 2. order genes: chromosome (natural), then start coordinate
##------------------------------------------------------------------
@gene_order = sort {
    natcmp($gene{$a}{f}[0], $gene{$b}{f}[0])
      || $gene{$a}{f}[3] <=> $gene{$b}{f}[3]
      || $gene{$a}{f}[4] <=> $gene{$b}{f}[4]
} @gene_order;

##------------------------------------------------------------------
## 3. assign new names
##------------------------------------------------------------------
open my $OUT, '>', $out     or die "$out: $!\n";
open my $MAP, '>', $mapfile or die "$mapfile: $!\n";
print $OUT "##gff-version 3\n";
print $MAP join("\t", qw(new_id old_id type chr start end strand)), "\n";

my (%counter, %seen_new);
my ($n_gene, $n_mrna) = (0, 0);

for my $gid (@gene_order) {
    my $g   = $gene{$gid};
    my @gf  = @{ $g->{f} };
    my $chr = $gf[0];

    # Chr5A -> chrom "5", hap "A"
    my ($num, $hap) = $chr =~ /^\D*(\d+)\s*([A-Za-z]*)$/;
    unless (defined $num) {
        warn "Cannot parse chromosome '$chr' - using it verbatim\n";
        $num = $chr; $hap = '';
    }
    my $tag = sprintf("%01d%s", $num, uc $hap);

    my $serial = ($counter{$tag} += $step);
    my $new_g  = sprintf("%s%sG%0*d", $prefix, $tag, $digits, $serial);
    die "Duplicate generated ID $new_g\n" if $seen_new{$new_g}++;

    $gf[1] = $source if defined $source;
    $gf[8] = "ID=$new_g;prev_id=$gid";
    $gf[8] .= ";Name=$1" if $keep_name && $g->{f}[8] =~ /\bName=([^;]+)/;
    print $OUT join("\t", @gf), "\n";
    print $MAP join("\t", $new_g, $gid, 'gene', @gf[0,3,4,6]), "\n";
    $n_gene++;

    # isoforms: longest CDS first, so .1 is the representative model
    my @iso = sort {
        cds_len($b) <=> cds_len($a)
          || $mrna{$a}{f}[3] <=> $mrna{$b}{f}[3]
          || $a cmp $b
    } @{ $g->{mrnas} };

    unless (@iso) { warn "gene $gid has no mRNA\n"; next }

    my $k = 0;
    for my $mid (@iso) {
        $k++;
        my $new_m = "$new_g.$k";
        my @mf = @{ $mrna{$mid}{f} };
        $mf[1] = $source if defined $source;
        $mf[8] = "ID=$new_m;Parent=$new_g;prev_id=$mid";
        print $OUT join("\t", @mf), "\n";
        print $MAP join("\t", $new_m, $mid, 'mRNA', @mf[0,3,4,6]), "\n";
        $n_mrna++;

        # children, in coordinate order; count each type separately
        my @kids = sort { $a->[3] <=> $b->[3] || $a->[4] <=> $b->[4] }
                   @{ $child{$mid} || [] };
        # for minus-strand genes exon1 should be the 5'-most exon
        @kids = reverse @kids if $mf[6] eq '-';

        my %n;
        for my $c (@kids) {
            my @cf = @$c;
            $cf[1] = $source if defined $source;
            my $t = $cf[2];
            my $sfx = $t eq 'exon'            ? 'exon'
                    : $t eq 'CDS'             ? 'cds'
                    : $t eq 'five_prime_UTR'  ? 'utr5p'
                    : $t eq 'three_prime_UTR' ? 'utr3p'
                    : lc $t;
            $cf[8] = sprintf("ID=%s.%s%d;Parent=%s", $new_m, $sfx, ++$n{$sfx}, $new_m);
            print $OUT join("\t", @cf), "\n";
        }
    }
}
close $OUT;
close $MAP;

warn "wrote $n_gene genes / $n_mrna mRNAs to $out\n";
warn "id map: $mapfile\n";

##------------------------------------------------------------------
sub cds_len {
    my ($mid) = @_;
    my $L = 0;
    for my $c (@{ $child{$mid} || [] }) {
        $L += $c->[4] - $c->[3] + 1 if $c->[2] eq 'CDS';
    }
    return $L;
}

sub natcmp {
    my @a = ($_[0] =~ /(\d+|\D+)/g);
    my @b = ($_[1] =~ /(\d+|\D+)/g);
    while (@a && @b) {
        my ($x, $y) = (shift @a, shift @b);
        my $c = ($x =~ /^\d+$/ && $y =~ /^\d+$/) ? $x <=> $y : $x cmp $y;
        return $c if $c;
    }
    return @a <=> @b;
}

sub usage {
    my ($m) = @_;
    print STDERR "ERROR: $m\n\n" if $m;
    print STDERR <<"END";
Usage:
  perl $0 -i CL_pasa.gff3 -o CL.renamed.gff3 -p Cla

  -i, --in     FILE   input GFF3 (.gz ok)
  -o, --out    FILE   output GFF3
  -m, --map    FILE   old->new table              [<out>.idmap.tsv]
  -p, --prefix STR    species prefix              [Cla]
      --step   INT    gene number increment       [10]
      --digits INT    digits in the serial        [5]
      --source STR    force column 2 to this value
      --keep-name     keep the original Name= attribute

Produces IDs like  Cla05AG00160  and  Cla05AG00160.1
END
    exit($m ? 1 : 0);
}
