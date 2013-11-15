#
# How to add new tests:
#
# Step 1:
# Create a new directory under the tests directory that matches the
# wildcard pattern: test*
#
# Step 2:
# In that new directory create 2 files, one file should match this 
# wildcard pattern: *input.txt
# and the other should match this wildcard pattern: *output.txt
# You can add other files to each directory so long as they do not
# match one of the wildcards, then they will be ignored.
#
# Step 3:
# In the input file add the commands you wish to run during the test
# These can be captured with the LOG command in NUI if you like.
#
# Step 4:
# In the output file put the contents of the output from NUI when the
# input commands are run.
#
# Step 5:
# Run make test from the main Morphy directory
#
# Step 6:
# If the output file does not exactly match the output NUI gave during
# the make test run, then the test failed and you need to fix it!
#
# To run the tests, you have the following params:
#   to update one specific test case:
#       UPDATE=test1
#   to update all test cases:
#       UPDATE=all
#   to run the "quick" tests
#       QUICKTEST=1
#
$(TESTDIRS): 
	$(eval TMP_IN := $(wildcard $@/*input.txt))
	$(eval TMP_OUT := $(wildcard $@/*output.txt))
	@./runtest.pl --test=$@ --input="$(TMP_IN)" --expected="$(TMP_OUT)" --update="$(UPDATE)" --quick="$(QUICKTEST)"


