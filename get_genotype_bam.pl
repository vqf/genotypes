#!/usr/bin/perl -w
use strict;

my $bam = shift;
my $chr = shift;
my $f = shift;
my $t = shift;
my $type = shift;
my $SAMTOOLS = shift || 'samtools';
if (!$bam || !$chr || !$f || !$t || !$type){
  die("Use: perl get_genotype_bam.pl bam_file chr from to type_file [samtools_executable=samtools]");
}

my $intypes = {};
open (IN, $type);
while (<IN>){
  chomp;
  my @f = split(/\s+/);
  my $el = shift @f;
  my $txt = shift(@f) || '0';
  $intypes->{$el} = $txt;
}
close IN;

my $header = `$SAMTOOLS view -H $bam`;

open (IN, "$SAMTOOLS view $bam '$chr:$f-$t' |");

my $out = '';
my $nlines = 0;
my $types = {};
my $topnlines = 100000;
while (<IN>){
  chomp;
  $nlines++;
  my @f = split(/\s+/);
  my $s = get_bam_seq(\@f, $f, $t);
  #print "$s\n";
  push @{$types->{$s}}, $_;
  last if ($topnlines && $nlines > $topnlines);
}
close IN;

my $counter = 0;
foreach my $t (keys %$types){
  if (exists($intypes->{$t})){
    my $fi = $intypes->{$t};
    if (!$fi){
      $fi = "bam_$counter.bam";
        $counter++;
    }
    my $nl = scalar @{$types->{$t}};
    warn("Writing to $fi, $nl lines\n");
    my $bamreads = join("\n", @{$types->{$t}});
    write_bam($fi, $header, $bamreads);
  }
}

sub write_bam{
  my $file = shift;
  my $h = shift;
  my $o = shift;
  open (OUT, "| $SAMTOOLS view -Shu - | $SAMTOOLS sort -o $file");
  print OUT $h . $o;
  close OUT;
  `$SAMTOOLS index $file`;
}

sub get_bam_seq{
  my $hit  = shift;
  my $from = shift;
  my $to   = shift;
  my $pos = $hit->[3];
  my $s   = $hit->[9];
  my $q   = $hit->[10];
  return '' if (length($s) < 10);
  my $result = '';
  my $opos = $pos;
  my $tpos = $pos;
  my $oseq = $s;
  my $pointer = 0;
  my $ins = 0;
  # Change sequence according to CIGAR (dels as D)
  my $cigar = $hit->[5];
  my @cig = $cigar =~ /(\d+\w)/g;
  return '' if ($cigar =~ /[H]/);
  my $cigpos = 0;
  foreach my $code (@cig){
    my ($ncode, $lcode) = $code =~ /(\d+)(\w)/;
    $lcode = "M" unless ($lcode =~ /[MIDNSHP]/);
    if ($lcode eq 'D'){
      my $tmp = '-' x $ncode;
      substr($s, $cigpos, 0, $tmp);
    }
    $cigpos += $ncode;
  }
  foreach my $code (@cig){
    my ($ncode, $lcode) = $code =~ /(\d+)(\w)/;
    $lcode = "M" unless ($lcode =~ /[MIDNSHP]/);
    if ($lcode eq 'M' || $lcode eq 'D' || $lcode eq 'N'){
      $tpos += $ncode;
      if ($tpos >= $from){
        if ($lcode eq 'D' || $lcode eq 'N'){
          $result = 'INDEL_DEL_'.$result;
        }
        my $start = max($pos, $from);
        my $end   = min($tpos, $to);
        my $ls = $start - $opos + $ins;
        my $ll = $end   - $start;
        $result .= substr($s, $ls, $ll);
      }
    }
    elsif ($lcode eq 'S'){  #Soft clipping
      $tpos += $ncode;
      $ins += $ncode;
      $result = 'S_' . $result;
      if ($tpos >= $from){
        return '';
      }
    }
    elsif ($lcode eq 'I'){  #Insertion
      $ins += $ncode;
      if ($tpos >= $from){
        $result = 'INDEL_INS_'.$result;
        my $toadd = substr($oseq, $pointer, $ncode);
        $result .= '(' . $toadd . ')'; 
      }
    }
    return $result if ($tpos >= $to);
    $pointer += $ncode unless ($lcode eq 'D' || $lcode eq 'N');
    $pos = $tpos;
  }
  return '';
}



sub min{
  my $result = shift;
  foreach my $n (@_){
    $result = $n if ($n < $result);
  }
  return $result;
}

sub max{
  my $result = shift;
  foreach my $n (@_){
    $result = $n if ($n > $result);
  }
  return $result;
}