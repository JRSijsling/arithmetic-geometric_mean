Description
-----------

This repository contains Pari and Sage code for calculating period matrices of genus 2 hyperelliptic curves by using the arithmetic-geometric mean. It will later be generalized to curves of higher genus.

Installation
------------

An installation of the computer algebra system Sage is needed to run this code.

It can be loaded inside Sage by going to the package/ directory and typing

load('Initialize.sage')

Alternatively, add the following lines to the Sage initialization file (typically found in ~/.sage/init.sage) to enable the functionality on startup anywhere in the system:

\_\_agmdir\_\_ = '~/[PATH]/package/'  
load(\_\_agmdir\_\_ + 'Initialize.sage')

where [PATH] is the path to the cloned or copied repository.

Usage
-----

An example that tests the routines in this package can be found in Example.sage. A debug file that interacts with the innards is in Debug.sage.
