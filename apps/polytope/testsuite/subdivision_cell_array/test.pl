my $c=load('1');
my @a=subdivision_cell_array($c);

compare_object('1a',$a[0])
and
compare_object('1b',$a[1])
and
    compare_object('1c',$a[2]);

