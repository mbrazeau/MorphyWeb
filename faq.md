---
title: FAQ
layout: default
filename: faq.md
---

Note: This version of the FAQ was written on 5 April 2020 for version 0.2 beta.

## Q: Why does the world need another phylogeny program?
Perhaps it doesn't. Morphy has its own way of encoding discrete data internally, it's own routines for deriving the score of a tree, and its own methods for accelerating the scoring of trees during searches. There is also a lack of serious parsimony programs available that are open-source and would have allowed me to invade the code and re-write these routines. There would be little hope I could have done this without my own practice at writing and building phylogeny programs. So now the world has Morphy. 

## Q: Why does Morphy only do parsimony?
Because likelihood calculations are slower and more complicated to program. I hope one day to get there.

## Q: How many states can Morphy handle?
Morphy can accommodate characters with up to 31 states plus the inapplicable token.

## Q: Does Morphy read TNT files?
At present, no. Please convert to a standard nexus file.

## Q: How can I contribute to Morphy?
There are several ways to contribute to this project. You can [submit bug reports](https://github.com/mbrazeau/morphy.archive/issues) when you encounter any issues. If you are a C or C++ developer, please get in touch about working on the project. 

## Q: Why can't Morphy do X?
Morphy is an immature program. There might be a plan to add that feature at some point in the future. If you'd like to see this feature added to Morphy, please contact the developers or [submit a feature request](https://github.com/mbrazeau/morphy.archive/issues).

## Q: Why is the Morphy 'business library' written in C while the user interface is written in C++?
The `mpl` library was written by Martin Brazeau who is most comfortable with the C language. The `NUI` interface for Morphy was mostly written by Chris Desjardins, who is a real programmer. He is smart and uses modern languages like C++. Probably one day the `mpl` will be re-written in C++.
