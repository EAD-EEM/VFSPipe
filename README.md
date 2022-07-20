# VFSPipe
**Tools for meshing PRZM, VFSMOD, and VVWM**

## Overview
The current version requires the user to set up a normal PWC2.0x run, and run it. Support for PWC1.52 is being phased out.

The user must provide the name of the .PWC file and the path to it and the outputs of that run (assumed to be the same for now). There are few direct user inputs, but these include a list of Strip Widths, a Shape Flag, a Ksat or choice of Ksat methods, and a choice of pesticide removal equations.

Executing the code runs VFSMOD for every runoff event, creates new .zts and VVWMTransfer files and runs VVWM. It also outputs a summary file of the inputs to and from the vegetative filter strip.

## Required External Files
The path to a folder containing the executables for PRZM, VFSMOD, and VVWM is required.
The list of soil texture classes is no longer required as it is coded into the project.

The path to a folder containing .zts files with the initial soil moisture conditions for the VFS for the scenario is still required.

## Theta .zts

VFSMod requires the initial soil moisture of the VFS at the start of each precipitation event, as well as the moisture content to calculate degradation between precipitation events. The very clever people at Waterborne figured out that you can get PRZM to output this in the .zts file. Since it is chemical independent, PRZM only needs to be run one time for each meteorological data file, so these .zts files are stored in a separate folder. Should scenarios be added, there's a julia function to make this easier.

## Future Improvements

Short Term
 - add checks for proper .zts file for theta (precedent moisture condition)
 - add output files comparing vvwm concentrations before and after VFSMOD run, both complete, and summary
 - add creation of the precedent moisture .zts into the main loop, since the PRZM run takes only a few seconds

Long Term
 - add ability to loop through scenarios, similar to a PWC scenario batch run
 - inclusion of the ability to write/edit PRZM files and run PRZM directly. This will allow automated, systematic exploration of more variables than are currently possible
 - a move to dynamic link library from text files to greatly speed execution. This requires modification of VFSMOD itself.
 - addition of a proper GUI
