# Morphy: Phylogenetic Analysis of Morphological Data With Character-state Inapplicability

## Download

The latest versions of Morphy [here](https://github.com/mbrazeau/morphy.archive/releases).

### Precompiled binaries:

[macOS](https://github.com/mbrazeau/morphy.archive/releases/download/0.2-beta/morphy_nui_v02b-macOS)

[Windows x64 (64-bit)](https://github.com/mbrazeau/morphy.archive/releases/download/0.2-beta/morphy_nui_v02b-win64.exe)

[Windows x86 (32-bit)](https://github.com/mbrazeau/morphy.archive/releases/download/0.2-beta/morphy_nui_v02b-win32.exe)

Linux users: follow the instructions [here](https://github.com/mbrazeau/morphy.archive) to clone all repos and run the build script.

## Install

On most systems, running the program shouldn't require any special installation. Mac users may need to follow the instructions below.

### Mac:
To install on macOS you will need to use the Terminal and navigate to the directory where the Morphy executable is stored and enter the command:

	sudo chmod a+x <morphy executable name>
	
To run Morphy, simply issue the command:

	./<morphy executable name>
	
Macs are quite restrictive about permissions for programs downloaded from the web. So, when you first try to run Morphy, the  operating system will try to stop you. Navigate to `System Preferences > Security & Privacy` and choose "Allow" for Morphy.

## Staying up to date

Join the [Morphy user group](https://groups.google.com/forum/#!forum/morphy-phylogenetic-software-announcements). This is a low-traffic forum for announcements about updates to the software and documentation.

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

## How do I cite Morphy?
There currently isn't a citable object for Morphy. Some of the speed-improving optimisations are still experimental and undergoing extensive testing. So, I am not sure I recommend using it for publication yet, but you are welcome to if you would like. For now, please cite [Brazeau, M.D., Guillerme, T. and Smith, M.R. 2019. An algorithm for Morphological Phylogenetic Analysis with Inapplicable Data. Systematic Biology. ](https://academic.oup.com/sysbio/article/68/4/619/5238046)
