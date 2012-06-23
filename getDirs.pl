#!/usr/bin/perl
use strict;
use Getopt::Long;

my $g_directory;

GetOptions( "directory=s" => \$g_directory);

sub getTestDirs(@)
{
    my ($directory) = @_;

    my @sorted= grep {s/(^|\D)0+(\d)/$1$2/g,1} sort grep {s/(\d+)/sprintf"%06.6d",$1/ge,1} <$directory/test*>;

    return @sorted;
}

my @dirs = getTestDirs($g_directory);
print("@dirs\n");

