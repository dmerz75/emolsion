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
* cscope - developer's tool for browsing program code. Run `./scope.sh`, creates cscope.files and cscope.out.
* [Boost] (http://www.boost.org/) - provides free peer-reviewed portable C++ source libraries.

## Microtubules' Contacts:
# The microtubule map: the diagram
            West
            4  7
         2  0  1  5
            3  6
            East

If 0 is alpha, 1 is beta,

then, 9 situations of 6 neighbors arise:

0    intra-alpha

1    intra-beta

0-1  intra-dimer

0-2  alpha-south

0-3  alpha-east

0-4  alpha-west

1-5  beta-north

1-6  beta-east

1-7  beta-west



## Notes/ Disclaimer/ Acknowledgements:
Six header files are not mine. They are currently used because I believe permission was given
provided the copyrights were retained and credit/notice given. They are:

* [Zed A. Shaw's] (https://learncodethehardway.org/c/) awesome debug macros.
  * debug.h
* These 5 headers are from UIUC, for use with their proprietary 'dcd' format. Information is
provided [here.] (http://www.ks.uiuc.edu/Research/vmd/plugins/doxygen/files.html)

        (C) Copyright 1995-2006 The Board of Trustees of the
                        University of Illinois
                         All Rights Reserved

  * dcdio.h ([dcdplugin.h](http://www.ks.uiuc.edu/Research/vmd/plugins/doxygen/dcdplugin_8c.html))
  * [endianswap.h] (http://www.ks.uiuc.edu/Research/vmd/plugins/doxygen/endianswap_8h.html)
  * [fastio.h] (http://www.ks.uiuc.edu/Research/vmd/plugins/doxygen/fastio_8h.html)
  * [largefiles.h] (http://www.ks.uiuc.edu/Research/vmd/plugins/doxygen/largefiles_8h.html)
  * [molfile_plugin.h] (http://www.ks.uiuc.edu/Research/vmd/plugins/doxygen/molfile__plugin_8h.html)

* Some of the coordinates in the test directory are straight from the [PDB] (http://www.rcsb.org/).
  * [2KHO.pdb] (http://www.rcsb.org/pdb/explore.do?structureId=2kho)
