#!/usr/bin/perl
use strict;
use Getopt::Long;
use Benchmark;

my $g_input_file;
my $g_expected_output_file;
my $g_test_dir;
my $g_update_dirs;
my $g_quick_test;

GetOptions( "input=s" => \$g_input_file,
            "expected:s" => \$g_expected_output_file,
            "test=s" => \$g_test_dir,
            "quick:s" => \$g_quick_test,
            "update:s" => \$g_update_dirs);

sub execute_nui(@)
{
    my ($input_file) = @_;
    open STDERR, '>&STDOUT';
    my $actual_results = `../../Debug/nui/nui < $input_file 2>&1`;
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
    my ($test_dir) = @_;
    my $tests_to_update;
    my @tests = split(/,/, $g_update_dirs);
    foreach my $test (@tests)
    {
        if (lc($test) eq "all")
        {
            $tests_to_update->{$test_dir} = 1;
            last
        }
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

sub print_results(@)
{
    my ($diffOutput, $diffSave, $dt) = @_;
    if ((length($diffOutput) == 0) && (length($diffSave) == 0))
    {
        printf("(dt=%-7.3f) - pass\n", $dt);
    }
    else
    {
        print("************** FAILURE\n$diffOutput\n\n$diffSave");
    }
}

sub process_results(@)
{
    my ($diffresults, $expected_output_file, $actual_output_file, $test_dir) = @_;
    if (length($diffresults) > 0)
    {
        my $tests_to_update = get_tests_to_update($test_dir);
        if ($tests_to_update->{$test_dir})
        {
            print("Updating test ");
            `cp $expected_output_file $expected_output_file.prev`;
            `cp $actual_output_file $expected_output_file`;
        }
    }
}

sub diff_results(@)
{
    my ($expected_output_file, $actual_output_file, $test_dir) = @_;
    my $diffresults = `diff $actual_output_file $expected_output_file 2>&1`;
    process_results($diffresults, $expected_output_file, $actual_output_file, $test_dir);
    return $diffresults;
}

sub runtest(@)
{
    my ($input_file, $expected_output_file, $test_dir, $quick_test) = @_;
    printf("%-8s", $test_dir);
    if (($quick_test) && (-e "$test_dir/longtest"))
    {
        print ("skip\n");
        return 0;
    }
    my $t0 = Benchmark->new;
    my $actual_results = execute_nui($input_file);
    my $t1 = Benchmark->new;
    my $actual_output_file = "$test_dir/actual.out";
    my $expected_save_file = "$test_dir/savefile_expected.txt";
    my $actual_save_file = "$test_dir/savefile.txt";

    $expected_output_file = update_expected_output($input_file, $expected_output_file);

    open(OUTPUTF, ">$actual_output_file");
    print OUTPUTF "$actual_results";
    close(OUTPUTF);
    my $dt = timediff($t1, $t0);
    my $diffOutput = diff_results($expected_output_file, $actual_output_file, $test_dir);
    my $diffSave = "";
    if (-e $expected_save_file)
    {
        $diffSave = diff_results($expected_save_file, $actual_save_file, $test_dir);
    }
    print_results($diffOutput, $diffSave, $dt->cpu_c);
    return 0;
}

exit(runtest($g_input_file, $g_expected_output_file, $g_test_dir, $g_quick_test));
