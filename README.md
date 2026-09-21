# PhotonArrivalFCS

MATLAB implementation of the Bayesian nonparametrics (BNP) method for analysing single-focus confocal data one photon at a time. The package collects the MATLAB scripts, auxiliary functions and the GUI implementing the algorithms described in the publication below.

## Reference

> Meysam Tavakoli, Sina Jazani, Ioannis Sgouralis, Omer M. Shafraz, Sanjeevi Sivasankar, Bryan Donaphon, Marcia Levitus and Steve Pressé,
> "Pitching Single-Focus Confocal Data Analysis One Photon at a Time with Bayesian Nonparametrics", *Physical Review X* **10**, 011021 (2020).
> <https://doi.org/10.1103/PhysRevX.10.011021>

Please cite the paper above in any publications using this software.

## Setup and Usage

1. To use this software in GUI form, run `PhotonArrivalFCS` in MATLAB's command window.
2. See [Help.pdf](Help.pdf) for an annotated walkthrough of the GUI.
3. For further details, or to report bugs in the original software, contact
   <mtavakol@purdue.edu> or <spresse@asu.edu>.

## MATLAB Versions

The source code and GUI were originally developed in MATLAB R2016a. They have since been updated and tested in MATLAB R2026a by Weiqing Xu.

The GUI was built with GUIDE, which has been removed from MATLAB, so the original code no longer ran on current releases. The changes rebuild the `handles` structure when the figure is opened on its own, and drop the references to controls that are no longer part of `PhotonArrivalFCS.fig`.

## Copyright and License

Copyright (C) 2019, Meysam Tavakoli, and Steve Presse

Permission is granted for anyone to copy, use, or modify these programs and accompanying documents for purposes of research or education, provided this copyright notice is retained, and note is made of any changes that have been made.

These programs and documents are distributed without any warranty, express or implied. As the programs were written for research purposes only, they have not been tested to the degree that would be advisable in any important application. All use of these programs is entirely at the user's own risk.
