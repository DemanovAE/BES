# Quick information
***

It is a simple example of scripts/macros sets suitable for work on the MEPhI (basov) cluster.

## Setting environment

Copy from a git repo:

        git clone https://devel.mephi.ru/PEParfenov/hpc_scripts.git

Change paths accordingly:

* 1 In `set_env.sh`: Change `ST_FEMTO_DST_INC_DIR` to standard one. It's the directory where `libStFemtoDst.so` is stored.
* 2 In `scripts/run.sh`: change path to a last line to one where your build will be.
* 3 In `scripts/start.sh`: change path to output directory (`OUTPUT_DIR`).

## Installation

Source required environment variables:

        cd hpc_scripts/
        . set_env.sh

Make new build directory. For example:

        mkdir build/
        cd build/

Generate makefile & install:

        cmake ../macro/
        make

Keep in mind that one has to do `make` command after changing `FemtoDstAnalyzer.C` in order to compile new code.

## Generate filelists

Use `GenerateLists.sh` to make filelists:

        . GenerateLists.sh FEMTODST_DIR N_FILES_IN_LIST

where `FEMTODST_DIR` - path to the directory with `femtoDst.root` files.
And `N_FILES_IN_LIST` denotes the maximum number of `femtoDst.root` files in each filelist.
Basic example:

          . GenerateLists.sh /mnt/pool/rhic/2/nigmatkulov/femtoDst/auau/200gev/12135/ 100

Resulting filelists will be in the `hpc_scripts/lists/` directory.

## Usage

### Interactive mode

To use `FemtoDstAnalyzer.C` in interactive mode:

        cd build/
        ./FemtoDstAnalyzer -i INPUTFILE -o OUTPUTFILE

where `INPUTFILE` - is input file or filelist with `femtoDst.root`.
`OUTPUTFILE` - is resulting root file.
Basic example:

        ./FemtoDstAnalyzer -i ../lists/StRuns1.list -o ./test.root

### Batch mode

To send jobs to basov cluster, use `scripts/start.sh`:

        . start.sh INPUT_FILELIST_DIR

where `INPUT_FILELIST_DIR` - is the directory where you store filelists.
Basic example:

        . start.sh /mnt/pool/rhic/4/parfenovpeter/STAR/Analysis/hpc_scripts/lists
# BES 
