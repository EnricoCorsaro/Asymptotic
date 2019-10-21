# Asymptotic extension to DIAMONDS

<p align="center">
<a href="https://github.com/EnricoCorsaro/Asymptotic"><img src="https://img.shields.io/badge/GitHub-Asymptotic-yellow"/></a>
<a href="https://github.com/EnricoCorsaro/Asymptotic/blob/master/LICENSE.txt"><img src="https://img.shields.io/badge/license-CC%20BY--SA-blue"/></a>
<a href='https://diamonds.readthedocs.io/en/latest/?badge=latest'><img src='https://readthedocs.org/projects/diamonds/badge/?version=latest' alt='Documentation Status' /></a>
<a href="https://github.com/EnricoCorsaro/Asymptotic/issues"><img src="https://img.shields.io/github/issues-closed/EnricoCorsaro/Asymptotic"/></a>
<img width="500" src="https://raw.githubusercontent.com/EnricoCorsaro/DIAMONDS/master/docs/figures/DIAMONDS_LOGO_WHITE.png"/>
</p>

### Author
- [Enrico Corsaro](mailto:enrico.corsaro@inaf.it)

### Short description
<div align="justify">
The Asymptotic extension to DIAMONDS is an extension of the Bayesian inference code DIAMONDS to perform the fitting of the asymptotic relations for p modes and for mixed modes by using as inputs a list of oscillation frequencies and corresponding uncertainties. The tool distinguished between radial modes and non radial modes, and between dipolar mixed modes and dipolar puressure modes. It operates in the directory containing the PeakBagging extension, because it assumes that the extracted
oscillation frequencies are first obtained through the PeakBagging code.
</div>

### Download & Installation
The procedure to retrieve the Asymptotic extension is identical to that of DIAMONDS (see [diamonds.readthedocs.io](http://diamonds.readthedocs.io/) for detailed information), so you can either clone the repository or simply download it as a ZIP file. In this second option, by unpacking the Asymptotic-master.zip file you will find a folder labeled Asymptotic-master and containing a structure similar to that of the folder Diamonds. First you need to rename the folder as **Asymptotic**, and place it in the same working directory of Diamonds (not inside the Diamonds folder!). This extension needs to be compiled separately from Diamonds, but only after you have compiled Diamonds first. Diamonds is used as a library for this extension. The compilation commands are the same as for Diamonds.  

**IMPORTANT**: Before proceeding with the compilation of the Asymptotic code make sure you put the Asymptotic folder at the same path level of that of Diamonds. This means that the Asymptotic folder has not to be placed inside the Diamonds folder, but inside the parent directory where you placed Diamonds.

### Documentation
Please make sure you read the documentation at [diamonds.readthedocs.io](http://diamonds.readthedocs.io/) before installing and using the code. This extension requires that the DIAMONDS code is first installed in your system. The installation of the Asymptotic extension is the same as that done for DIAMONDS.

### Tutorials
To run the tutorials provided in the package, please follow the guidelines presented in [tutorials/README.md](https://github.com/EnricoCorsaro/Asymptotic/blob/master/tutorials/README.md)
