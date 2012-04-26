#!/usr/bin/perl
use strict;
use Getopt::Long;

my $g_input_file;
my $g_expected_output_file;
my $g_test_dir;
my $g_update_dirs;
GetOptions( "input=s" => \$g_input_file,
            "expected:s" => \$g_expected_output_file,
            "test=s" => \$g_test_dir,
            "update:s" => \$g_update_dirs);

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

sub get_tests_to_update(@)
{
    my $tests_to_update;
    my @tests = split(/,/, $g_update_dirs);
    foreach my $test (@tests)
    {
        $tests_to_update->{$test} = 1;
    }
    return $tests_to_update;
}

sub update_expected_output(@)
{
    my ($input_file, $expected_output_file) = @_;
    if (length($expected_output_file) == 0)
    {
        $expected_output_file = $input_file;
        $expected_output_file =~ s/input\.txt/output\.txt/g;
    }
    return $expected_output_file;
}

sub process_results(@)
{
    my ($diffresults, $input_file, $expected_output_file, $actual_output_file, $test_dir) = @_;
    if (length($diffresults) == 0)
    {
        print ("  $test_dir - pass\n");
    }
    else
    {
        my $tests_to_update = get_tests_to_update();
        print("\n\n\n************* $test_dir FAILURE - $diffresults\n");
        if ($tests_to_update->{$test_dir})
        {
            print("Updating test: $test_dir\n\n");
            `cp $expected_output_file $expected_output_file.prev`;
            `cp $actual_output_file $expected_output_file`;
        }
    }
}

sub runtest(@)
{
    my ($input_file, $expected_output_file, $test_dir) = @_;
    my $actual_results = execute_nui($input_file);
    my $actual_output_file = "$test_dir/actual.out";

    $expected_output_file = update_expected_output($input_file, $expected_output_file);

    open(OUTPUTF, ">$actual_output_file");
    print OUTPUTF "$actual_results";
    close(OUTPUTF);
    my $diffresults = `diff $actual_output_file $expected_output_file 2>&1`;

    process_results($diffresults, $input_file, $expected_output_file, $actual_output_file, $test_dir);

    return 0;
}

exit(runtest($g_input_file, $g_expected_output_file, $g_test_dir));