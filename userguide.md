## User guide

A list of commands is provided with the ‘help’ command. However, that list is a lot of promises and not a lot of action. The basic functionality of Morphy is doing searches using either random addition sequence + branch-swapping or using the "ratchet”.

Basic command structure:

	<command> = <option>

### Changing your working directory:

	cd=/path/to/working/directory

Opening a file:

	open=mynexusfile.nex

Running a simple heuristic search using the default settings and input order of the taxa:

	heuristic

or simply ‘heur’ is enough to initiate a search.

We know that’s not a very good approach, so we’ll want to set some additional parameters, do multiple random addition sequences and branch-swappers.

The following commands are global options and persistent.

Set the addition sequence:

	addseq=random (other options: asis [the default]; none [swaps trees in memory])

Set the number of replicates:

	nreps=N (where N is he number of reps you wish to run)

Then run the ‘heur’ command and your RAS+TBR search will begin. (You can also set branchswap=spr if you want to do an SPR search which is a bit faster, but less thorough).

This will be more effective at finding optimal trees. However, large datasets with lots of inapplicable data in them will tend to run more slowly as the program needs to run a larger number of full downpasses to get the length of the tree (most characters are checked locally which is very fast, but when that’s not possible a full downpass needs doing). This slows the program quite a bit making it less efficient.

Another issue is the problem of optimality islands: for whatever reason, this algorithm defines more suboptimal islands than standard parsimony and is therefore more prone to getting ’trapped’ on a local optimum. Combine that with overall lower efficiency of the length-counting algorithm and you can have a very slow search! Luckily, the “ratchet” is extremely effective at breaking out of these local optima and allowing you to run more replicates faster. 

To run a ratchet search, you simply need to use the command:

	ratchet=yes

It is highly recommended that you are using TBR for this, however.

After the ratchet search completes, you can do full TBR or SPR search on the trees in memory. So, an example workflow might be as follows:

	open=mydata.nex
	addseq=random
	nreps=100
	ratchet=yes
	heur

… wait …

	addseq=none
	ratchet=no
	heur

… wait even longer …

The search will complete eventually (or you can hit command+period at any time to break out of the search).

At present, Morphy doesn’t compute a consensus tree or output a text-drawing of your tree(s). You’ll want to export your trees to a tree file to view them or compute the consensus.
That’s simply done with:

	save=mytrees.tre 

Or whatever you want to call your file, of course. Note: this currently overwrites without warning.

