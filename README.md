# Emolsion
A pdb and dcd reading/writing molecular mechanics evaluation program.

## Quickstart:
1. DCD reading: To read a pdb with corresponding dcd (from frames 6 - 27, stepping by 3 frames) ..

        make dcdr
        emol_dcdreader mt.ref.pdb mt_partial.dcd 6 27 3

2. DCD writing: To write a dcd from the pdb with corresponding dcd (from frames 4 - 30,
stepping by 2 frames)..

        make dcdw
        emol_dcdwriter mt.ref.pdb mt_partial.dcd 4 30 2



## Dependencies:
Other softwares I often use with this project include:
* cscope - developers tool for browsing program code
* [Boost] (http://www.boost.org/) - provides free peer-reviewed portable C++ source libraries.



## Disclaimer:
Six header files are not mine. They are currently used because I believe permission was given
provided the copyrights were retained and credit/notice given. They are:

* [debug.h] (https://learncodethehardway.org/c/) is from Zed A. Shaw.
* These 5 headers are from UIUC, for use with their proprietary 'dcd' format.

        (C) Copyright 1995-2006 The Board of Trustees of the
                        University of Illinois
                         All Rights Reserved

  * dcdio.h
  * endianswap.h
  * fastio.h
  * largefiles.h
  * molfile_plugin.h

Some of the coordinates in the test directory are straight from the [PDB] (http://www.rcsb.org/).
* [2KHO.pdb] (http://www.rcsb.org/pdb/explore.do?structureId=2kho)
