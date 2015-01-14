#! /usr/bin/perl -w

my $method= 'camB3LYP';

my $molecule =
'molecule_inline
C    0.0000    0.0000    0.0000
O    2.0   0.0   0.0
O   -2.0   0.1   0.0
EOF
basis="3-21G"';

sub runField($) {
    my $field =$_[0];
    my $res = 0;
    unlink "ergoscf.out";
    open ERGO, "|source/ergo" or die "Cannot exec ergo: $!\n";
    print ERGO << "EOI";
$molecule
scf.electric_field_z = $field
run "$method"
EOI
    close ERGO or die "ergo execution failed: $!\n";
    open OUT, "ergoscf.out" or die "Cannot open ergoscf.out for reading: $!\n";
    while(<OUT>) {
        if(/FINAL/){
            print;
            $res = (split ' ', $_)[3];
        }
    }
    close OUT or die "Failed to close file: $!\n";
    return $res;
}

sub runPola {
    unlink "ergoscf.out";
    open ERGO, "|source/ergo" or die "Cannot exec ergo: $!\n";
    print ERGO << "EOI";
$molecule
get_polarisability "$method" "Z" 0.0
EOI
    close ERGO or die "ergo execution failed: $!\n";
    open OUT, "ergoscf.out" or die "Cannot open ergoscf.out for reading: $!\n";
    while(<OUT>) {
        if(/^REMA Response << Z/) {
            print;
            $res = (split ' ', $_)[10];
        }
    }
    close OUT or die "Failed to close file: $!\n";
    return $res;
}

my $field = 1e-4;
my $eFieldZero  = runField 0;
my $eFieldMinus = runField -$field;
my $eFieldPlus  = runField +$field;
my $pol         = runPola;

print "Finite ", -(2*$eFieldZero-$eFieldPlus-$eFieldMinus)/($field*$field),
    " analytic: $pol\n";
