# This takes the .zts output from PRZM5 to parameterize the input files of VFSMOD
#for each individual precipitation event, and compiles the results into
#another .zts file upon which VVWM can be run
include("./VFSUtils.jl")
#include("C:\\Users\\hcuser\\github\\PRZMTools.jl\\src\\PRZMTools.jl")
#Can save some work by following a rigid filename convention
projectName = "SOILS" #Six characters only
stripWidthInM = 10.0
#Can replace the last two letters of the input files for the temporary working files so they do not get written over
twoCharacterCode = "RG"
workingPath = "Z:\\SharedwithVM\\VFS\\Soils\\9005\\"
pwcName = "Persistane-REG"
# If there's a proper hourly precipitation file, use that, otherwise generate the precipitation event programmatically from the daily value
useHPF = false
stormLengthInHours = 8 #only read when useHPF is false

#Less often changed are the paths to unchanging files
exePath = "Z:\\SharedwithVM\\VFS\\executables\\"
turfPath = "Z:\\SharedwithVM\\VFS\\CanadianTurfZts\\"

# Use PRZMTools to set up and run PRZM and VVWM
#S = PRZMTools.
#P = PRZMTools.Pesticide(string(pwcName, ".PWC")) # Pesticide from PWC
#C = PRZMTools.readCrop("Potato","Atlantic",datafile="C:\\Users\\hcuser\\github\\PRZMTools.jl\\Data\\cropdata.txt")

usInp = userInputs(projectName = projectName, stripWidthInM = stripWidthInM, twoCharacterCode = twoCharacterCode, workingPath = workingPath, pwcName = pwcName, useHPF = useHPF, stormLengthInHours = stormLengthInHours, exePath = exePath, turfPath = turfPath)

#This runs VFSMOD for all precipitation events, rewrites the zts file and returns the string to run VVWM
newVVWM, oldTXT, newTXT = vfsMain(usInp)
run(newVVWM)

noVFS = getVVWMText(columnNames(),oldTXT)
yesVFS = getVVWMText(columnNames(),newTXT)
plot(noVFS.Run, [noVFS.Peak,yesVFS.Peak], title = "Efficacy of VFS for EEC Reduction", label = ["Without VFS" "With VFS"], linecolor = ["Brown" "Green"],lw = 2, xlabel = "Simulation Year", ylabel = "PPB")
plot(noVFS.Run, [noVFS.OneYear,yesVFS.OneYear], title = "Efficacy of VFS for EEC Reduction", label = ["Without VFS" "With VFS"], linecolor = ["Brown" "Green"],lw = 2, xlabel = "Simulation Year", ylabel = "PPB")
