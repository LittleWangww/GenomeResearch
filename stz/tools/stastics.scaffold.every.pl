#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
my $fasta_file=shift;
my $opfile=shift;
my $catchseq_obj=Bio::SeqIO->new(-file=>$fasta_file,-format=>'fasta');
my %totalbase;
my @seqArray;
my $seq_N=0;
my $lengthall=0;
while(my $seq=$catchseq_obj->next_seq){
    my $baseseq=$seq->seq;
    chomp $baseseq;
    push(@seqArray,length($baseseq));
    $lengthall+=length($baseseq);
    $seq_N++;
    my @base=split(//,$baseseq);
    $totalbase{$_}++ foreach @base;
    }
}
my @seqArrayed=sort{$b<=>$a}@seqArray;
my $length50=$lengthall*(0.5);
my $length60=$lengthall*(0.6);
my $length70=$lengthall*(0.7);
my $length80=$lengthall*(0.8);
my $length90=$lengthall*(0.9);
my $length5=0;
my $N50;
foreach $_(@seqArrayed){
    $length5+=$_;
    if($length5>=$length50){
	$N50=$_;
	last;
    }
}
my $length6=0;
my $N60;
foreach $_(@seqArrayed){
    $length6+=$_;
    if($length6>=$length60){
	$N60=$_;
	last;
    }
}
my $length7=0;
my $N70;
foreach $_(@seqArrayed){
    $length7+=$_;
    if($length7>=$length70){
	$N70=$_;
	last;
    }
}
my $length8=0;
my $N80;
foreach $_(@seqArrayed){
    $length8+=$_;
    if($length8>=$length80){
	$N80=$_;
	last;
    }
}
my $length9=0;
my $N90;
foreach $_(@seqArrayed){
    $length9+=$_;
    if($length9>=$length90){
	$N90=$_;
	last;
    }
}
my %hash_N50;
my $n5=0;
my %hash_N60;
my $n6=0;
my %hash_N70;
my $n7=0;
my %hash_N80;
my $n8=0;
my %hash_N90;
my $n9=0;
while(my $seq2=$catchseq_obj->next_seq){
    my $baseseq2=$seq2->seq; chomp $baseseq2;
my @base2=split(//,$baseseq2);
if(length($baseseq2)==$N50){
    $n5++;
    $hash_N50{$_}++ foreach @base2;

}
elsif(length($baseseq2)==$N60){
    $n6++;
    $hash_N60{$_}++ foreach @base2;

}
elsif(length($baseseq2)==$N70){
    $n7++;
    $hash_N70{$_}++ foreach @base2;

}
elsif(length($baseseq2)==$N80){
    $n8++;
    $hash_N80{$_}++ foreach @base2;

}
elsif(length($baseseq2)==$N90){
    $n9++;
    $hash_N90{$_}++ foreach @base2;

}
}
my $n50_N=(($hash_N50{N}+$hash_N50{n})/$n5);
my $per_n50_N=($n50_N/$lengthall)*100;
my $n60_N=(($hash_N60{N}+$hash_N60{n})/$n6);
my $per_n60_N=($n60_N/$lengthall)*100;
my $n70_N=(($hash_N70{N}+$hash_N70{n})/$n7);
my $per_n70_N=($n70_N/$lengthall)*100;
my $n80_N=(($hash_N80{N}+$hash_N80{n})/$n8);
my $per_n80_N=($n80_N/$lengthall)*100;
my $n90_N=(($hash_N90{N}+$hash_N90{n})/$n9);
my $per_n90_N=($n90_N/$lengthall)*100;
my $AVG_length=int($lengthall/$seq_N);
my $Min_length=shift(@seqArray);
my $Max_length=pop(@seqArray);
my $GC=(($totalbase{G}+$totalbase{C})/$lengthall)*100;
open(O,">$opfile");
print O "Total_num\tTotal_length\tAVG_length\tMin_length\tMax_length\tGC%\n";
print O "$seq_N\t$lengthall\t$AVG_length\t$Min_length\t$Max_length\t$GC\n";
print O "N50\tN50_N\tN50_N%\tN60\tN60_N\tN60_N%\tN70\tN70_N\tN70_N%\tN80\tN80_N\tN80_N%\tN90\tN90_N\tN90_N%\t";
print O "$N50\t$n50_N\t$per_n50_N\t$N60\t$n60_N\t$per_n60_N\t$N70\t$n70_N\t$per_n70_N\t$N80\t$n80_N\t$per_n80_N\t$N90\t$n90_N\t$per_n90_N\t";
close O;
