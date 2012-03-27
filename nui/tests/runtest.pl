#!/usr/bin/perl
use strict;
use Getopt::Long;

my $g_input_file;
my $g_expected_output_file;
my $g_test_dir;
GetOptions( "input=s" => \$g_input_file,
            "expected=s" => \$g_expected_output_file,
            "test=s" => \$g_test_dir);

sub execute_nui(@)
{
    my ($input_file) = @_;
    my $actual_results = `../obj/nui < $input_file 2>&1`;
    return $actual_results;
}

sub read_expected_results(@)
{
    my ($expected_output_file) = @_;
    open FILE, $expected_output_file or die "Couldn't open file: $expected_output_file"; 
    my $expected_results = join("", <FILE>); 
    close FILE;
}

sub runtest(@)
{
    my ($input_file, $expected_output_file, $test_dir) = @_;
    my $actual_results = execute_nui($input_file);

    open(OUTPUTF, ">$test_dir/actual.out");
    print OUTPUTF "$actual_results";
    close(OUTPUTF);
    my $diffresults = `diff $test_dir/actual.out $expected_output_file`;
    if (length($diffresults) == 0)
    {
        print ("$test_dir - pass\n");
    }
    else
    {
        print("\n\n\n************* $test_dir FAILURE - $diffresults\n");
    }
    return 0;
}

exit(runtest($g_input_file, $g_expected_output_file, $g_test_dir));
