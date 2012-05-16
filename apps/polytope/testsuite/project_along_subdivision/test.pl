prefer_now "cdd";

my $c=load('c');

compare_object('1',project_along_sequence($c,new Array<Int>([1,1,3,3])));
