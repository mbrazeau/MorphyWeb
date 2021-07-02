---
title: Morphy Phylogenetic Library (mpl)
layout: default
filename: fxnlib.md
---
This is the documentation for the Morphy API. 
If you're not writing programs for phylogenetic analysis, then this page is not likely to be of interest to you.
If you're looking for MorphyLib, a C library that can be used for tree score calculations, download and documentation are available [here](https://github.com/mbrazeau/morphylib).
MorphyLib and the MPL are not meant to be used together in the same program.

**Introduction**

To facilitate using Morphy from external programs and environments, we have written Morphy as a library (MPL: Morphy Phylogenetic Library). 
It is therefore possible to build Morphy into your own programs with minimal effort and to also call its functions from environments such as Python or R. 


**Building the MPL**
<void>

**Using the MPL**
The API is [documented]() in the [`mpl.h`]() file. 

The caller creates an instance of the Morphy object called the `mpl_handle` through the `mpl_handle_new` function. No direct access to the members of this structure are required, and all interaction is through functions declared in `mpl.h`.

```C
#include "mpl.h"

...
mpl_handle handl;
handl = mpl_handle_new();
```

Morphy has its own rules for 'packing' state data as set bits in integers and therefore you are required to supply a data matrix and dimensions to it in order to perform any calculations.


```C
int ntax = 5;
int nchar = 10;

mpl_set_dimensions(ntax, nchar, handl); // Sets the dimensions of the matrix.
```
