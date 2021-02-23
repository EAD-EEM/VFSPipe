# This takes the .zts output from PRZM5 to parameterize the input files of VFSMOD
#for each individual precipitation event, and compiles the results into
#another .zts file upon which VVWM can be run
include("./VFSUtils.jl")
#include("C:\\Users\\hcuser\\github\\PRZMTools.jl\\src\\PRZMTools.jl")
#User choice
projectName = "ERIC" #Six characters only
stripWidthInM = 10.0
#User choice - Can replace the last two letters of the input files for the temporary VFSMOD working files so they do not get written over
#A holdover from VFSMOD which may have other applications
twoCharacterCode = "QQ"

# The path to the working directory - should contain, from PWC: .SWI or .PWC PWC file, .SCN scenario file, .dvf or .hpf weather file, and the output files from the original run (.zts, .inp, vvwmtransfer.txt)
workingPath = "Z:\\SharedwithVM\\VFS\\ExampleCase\\"
# Name of the .SWI or .PWC file
pwcName = "GenericChemical"
# If there's a proper hourly precipitation file, use that, otherwise generate the precipitation event programmatically from the daily value
useHPF = false
stormLengthInHours = 8 #only read when useHPF is false

#Less often changed are the paths to unchanging files
# Executables file must contain vfsm,exe and vvwm.exe, as well as SoilTriangle.csv
exePath = "Z:\\SharedwithVM\\VFS\\executables\\"

# The location of the pre-run .zts files where θ, the  pre-storm water content of the VFS, is found.
turfPath = "Z:\\SharedwithVM\\VFS\\TurfZts\\"
pesticideEquation = 3 # 1:Sabbagh;2:refitSabbagh;3:mass-bal.;4:Chen

# For exploration only - set to any negative value to use the lookup table method (based on soil texture pulled from the PWC scenario)
Ksat = -9999 #in m/s
    
# The dimensions of the VFS are calculated from its width (which VFS calls its length) and the area of the pond into which it flows
# It can be a round pond with a round VFS around it, a square pond with a square VFS around it, or a square pond with a rectangular VFS on one side
shapeFlag = 2 # 1 for round, 2 for square, 3 for rectangular

usInp = userInputs(projectName = projectName, stripWidthInM = stripWidthInM, twoCharacterCode = twoCharacterCode, workingPath = workingPath, pwcName = pwcName, useHPF = useHPF, stormLengthInHours = stormLengthInHours, exePath = exePath, turfPath = turfPath, pesticideEquation = pesticideEquation, Ksat=Ksat, shapeFlag = shapeFlag)

#This runs VFSMOD for all precipitation events, rewrites the zts file and returns the string to run VVWM
newVVWM, oldTXT, newTXT = vfsMain(usInp)
run(newVVWM)
noVFS = getVVWMText(columnNames(),oldTXT)
yesVFS = getVVWMText(columnNames(),newTXT)
plot(noVFS.Run, [noVFS.Peak,yesVFS.Peak], title = "Efficacy of VFS for EEC Reduction", label = ["Without VFS" "With VFS"], linecolor = ["Brown" "Green"],lw = 2, xlabel = "Simulation Year", ylabel = "PPB")
plot(noVFS.Run, [noVFS.OneYear,yesVFS.OneYear], title = "Efficacy of VFS for EEC Reduction", label = ["Without VFS" "With VFS"], linecolor = ["Brown" "Green"],lw = 2, xlabel = "Simulation Year", ylabel = "PPB")
