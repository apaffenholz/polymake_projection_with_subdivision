enable_if_configured("javaview.rules") or return;
script("test_filters");

prefer_now 'cdd';
prefer_now 'normaliz2';
my $c=load("c");

javaview($c->VISUAL->SUBDIVISION, File=>diff_with("1", filter_JVX()));
