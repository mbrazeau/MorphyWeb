---
title: User guide
layout: default
filename: userguide.md
---

## User guide

## Table of contents

Coming soon...

## Commands

A list of commands is provided with the ‘help’ command. However, that list is a lot of promises and not a lot of action. The basic functionality of Morphy is doing searches using either random addition sequence + branch-swapping or using the "ratchet”.

Basic command structure:

	<command> = <option>

## Opening files, setting working directory 

To set the working directory, the command is issued as follows:

	cd=/path/to/working/directory

Opening a file:

	open=mynexusfile.nex

## Input format

Morphy currently only reads Nexus files.

The following are valid state symbols:

	0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz

The following symbols can be used as 'wildcards' 

* `-` The gap symbol (can be any of inapplicable, missing, or an additional state). The options are global for the data set and are set using the `gap=` command, where the options are `inapplicable` (default), `missing`, and `state`.
* `?` Signifies missing data in all categories (applicable or inapplicable)
* `+` Signifies 'unknown' defined here to be: "applicable, but of unknown state". This could be, for instance, the colour, of a structure known in a fossil.

## Handling the gap state:

By default, Morphy treats the gap symbol in the data as inapplicable under the algorithm by Brazeau, Guillerme, and Smith (2019). However, it is possible to set the gap symbol to be treated either as missing or as an extra state. Simply use the `gap` command with either `inapplicable`, `missing`, or `state`.

## Setting the parsimony type:

Characters are unordered by default. To set characters to ordered parsimony, simply use the `ordered` command followed by an equals sign and the space-separated list of character numbers (the first character 1). For instance:

    ordered=1 7 21 30

To set these back, simply issue the inverse unordered command:

    unordered=1 7 21 30

At present, inputting ranges is not supported. A future version will allow setting ranges as well as using an 'all' declaration.

## Searching
Running a simple heuristic search using the default settings and input order of the taxa:

	heuristic

or simply ‘heur’ is enough to initiate a search.

We know that’s not a very good approach, so we’ll want to set some additional parameters, do multiple random addition sequences and branch-swappers.

The following commands are global options and persistent.

Set the addition sequence:

	addseq=random (other options: asis [the default]; none [swaps trees in memory])

It is possible to employ a 'hold' criterion, as in PAUP, which holds the top `N` most parsimonious trees at each step in the addition sequence. 

	hold=N

The default hold is 1. Setting higher hold values can prevent the frequency with which Morphy is entrapped by a local optimum.

Set the number of replicates:

	nreps=N (where N is he number of reps you wish to run)

Then run the ‘heur’ command and your RAS+TBR search will begin. (You can also set branchswap=spr if you want to do an SPR search which is a bit faster, but less thorough).

TBR will be more effective at finding optimal trees. 
However, large datasets with lots of inapplicable data in them will tend to run more slowly as the program needs to run a larger number of more complete downpasses to get the length of the tree (most characters are checked locally which is very fast, but when that’s not possible provisional updates have to be conducted over larger parts of the tree). 
This slows the program quite substantially.

Another issue is the problem of optimality islands: for whatever reason, this algorithm defines more suboptimal islands than standard parsimony and is therefore more prone to getting ’trapped’ on a local optimum. Combine that with overall lower efficiency of the length-counting algorithm and you can have a very slow search! Luckily, the “ratchet” is extremely effective at breaking out of these local optima and allowing you to find optimal trees much faster.

To turn on ratchet searching, simply issue the following command before searching:

	ratchet=yes

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

The search will complete eventually (You can interrupt a search with Cmd+period on a Mac or Ctrl+C on other systems).

## Saving results
At present, Morphy doesn’t compute a consensus tree or output a text-drawing of your tree(s). You’ll want to export your trees to a tree file to view them or compute the consensus.
That’s simply done with:

	save=mytrees.tre 

Or whatever you want to call your file, of course. 

