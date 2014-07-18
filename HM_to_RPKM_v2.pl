use strict;

### This script uses a HashMatch output file and a file containing transcript sequences to generate RPKM values

# ARGV[0]
# HashMatch output file
# Sb09g029520.1   1479    1514    -       GAGAATGCGATTTCCTCTAATGTTCTTCCCTTGGTC    >4732268        >6065878
# hash the read count info
my %count = ();
my %mapped = ();
open (FILE1, "<$ARGV[0]") || die "$!: exiting\n";
while (<FILE1>) {
    my @genes;
    my @data = split(/\t/,$_);
    my $transcript = shift(@data);
    my $start = shift(@data);
    my $stop = shift(@data);
    my $strand = shift(@data);
    my $sequence = shift(@data);
    # apply DUST to eliminate junk reads with homopolymer runs [kmay be too stringent]
    $sequence =~ s/((.+)\2{4,})/'N' x length $1/eg;
    if ($sequence =~ /N/) {
    }
    else {
	foreach my $item (@data) {
	    push @genes, $transcript;
	    $mapped{$item} = 1;
	}

	foreach my $gene (@genes) {
	    $count{$gene}++;
	}
    }
}
close FILE1;

# hash the transcript lengths
my $count = 0;
my ($one,$two,$three,$function,$homology,$seq) = "";
my %kbs = ();
open (FILE2, "<$ARGV[1]") || die "$!: exiting\n";
while (<FILE2>) {
    chomp($_);
    $count++;
    if ($_ =~ /^>/) {
        if ($count > 1) {
            my $kb = length($seq)/1000;
            $kbs{$one} = $kb;
            $one = $two = $three = $function = $homology = $seq = "";
        }
        my ($three_ids, $function, $homology) = split(/\t/, $_);
        $three_ids =~ s/>//;
        ($one, $two, $three) = split(/\#/, $three_ids);
        if (!$three) {
            $three = "NULL";
        }
    }
    else {
        $seq = $seq.$_;
    }
}
my $kb = length($seq)/1000;
$kbs{$one} = $kb;
close FILE2;

# get the total mapped read count;
my $mapped = keys(%mapped)/1000000;
# print "$mapped\n";

foreach my $gene (
sort { $count{$b} <=> $count{$a} } keys %kbs
  ) {
    my $length = $kbs{$gene};
    my $rpkm = ($count{$gene}/$length)/$mapped;
    
    print "$gene\t$rpkm\n";
  }
