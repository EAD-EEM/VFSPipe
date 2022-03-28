# VFSPipe
**Tools for meshing PRZM, VFSMOD, and VVWM**

## Overview
The current version requires the user to set up a normal PWC run, and run it.
The user must provide the name of the .SWI file and the path to it and the outputs of that run (assumed to be the same for now). There are few direct user inputs, but these include a list of Strip Widths, a Shape Flag, a Ksat or choice of Ksat methods, and a choice of pesticide removal equations.

Executing the code runs VFSMOD for every runoff event, creates new .zts and VVWMTransfer files and runs VVWM (and displays a graph if run in the REPL).

## Required External Files
The path to a folder containing the executables for VFSMOD and VVWM is required. (WOMM)
The list of soil texture classes is also required. It is not a large file and perhaps should be moved into the project.

The path to a folder containing .zts files with the initial soil moisture conditions for the VFS for every scenario is required. (WOMM)

## Turf .zts

VFSMod requires the initial soil moisture of the VFS at the start of each precipitation event. The very clever people at Waterborne figured out that you can get PRZM to output this in the .zts file. Since it is chemical independent, PRZM only needs to be run one time for each scenario, so these .zts files are stored in a separate folder. Should scenarios be added, there's a julia function to make this easy.

## Future Improvements

Short Term
 - add slope(s) as an optional user input (currently slope is read from the PWC scenario).

Long Term
 - inclusion of the ability to write/edit PRZM files and run PRZM directly. This will allow automated, systematic exploration of more variables than are currently possible
 - a move to dynamic link library from text files to greatly speed execution. This requires modification of VFSMOD itself.
 - addition of a proper GUI
