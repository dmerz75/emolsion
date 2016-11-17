# Emolsion
A pdb and dcd reading/writing molecular mechanics evaluation program.

## Quickstart:
To run and see an example:

    make dcd1

To run later (the executables currently write to the test directory):

    run_readpdb_def mt.ref.pdb mt_partial.dcd 6 27 3

where 6, 27, 3 are "start, stop, step," respectively, for the frames in a dcd (trajectory) file.






## Dependencies:
Other softwares I often use with this project include:
* cscope (extra/cscope 15.8.a-3) A developers tool for browsing program code
* [Boost] (http://www.boost.org/)

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
