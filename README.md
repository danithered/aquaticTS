# Evolution of thermal reaction curves in aquatic ecosystems

## Description

Simulation developed in **C++**, visualisations in **R**.

blahbla

To get a wider description about the model, or the program itself read the article (...) or visit the code documentation. For the latter, please move to the repository (`cd path-to-repo`) and run `doxygen`
To read the documentation open `path-to-repo/doc/html/index.html`

## Prequisites

Developed and tested on Debian 12.

Before installation, please make sure you have the followings installed and in `$PATH`.

- gcc, pkg-config, make (`sudo apt install build-essential`)
- GNU Scientific Library (https://www.gnu.org/software/gsl/)
- Boost (`sudo apt install libboost-all-dev`)

And of course you have to download the neccessary files. The easiest way to do so is:

1. `cd mypath`, where *mypath* is the directory where you wish to download the files to
2. `git clone git@github.com:danithered/aquaticTS.git`
3. the files will be in *mypath/aquaticTS/*

## Installation

1. `cd` to project directory
2. `make`

## Usage

### Running simulations

To just try it: the program completes with default parameters, so just run: `./simulation`. It will output into the `./OUT/test/` directory.

To get help with setting the paramters of the simulations see: `./simulation --help`

### Download climate data with the help of **R**

The file `./src/get_climate.R` gives you all the help you need to get input climate files for the simulations. As a prequisite you have to install **R** (`sudo apt install r-base`) and get a code editor for *.R* files (for biologists *RStudio* is recommended).

## Citation

soon :D
