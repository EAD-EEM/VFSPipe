# This takes the .zts output from PRZM5 to parameterize the input files of VFSMOD
# for each individual precipitation event, and compiles the results into
# another .zts file upon which VVWM is run
include("./VFSUtils.jl")
#include("C:\\Users\\hcuser\\github\\PRZMTools.jl\\src\\PRZMTools.jl")

#User choice
stripWidthInM = [11.0, 5.0] #[50.0, 15.0, 10.0, 5.0, 3.0]

# The path to the working directory - should contain the .SWI or .PWC PWC file, and the output files from the original run (.zts, .inp, vvwmtransfer.txt).
# If an .SCN file using the filename specified in the PWC run is located here, it will supersede one in the folder specified by scenarioPath
workingPath = "Z:\\SharedwithVM\\VFS\\Phase Two\\PipetestZ\\" # must end in a double back-slash or slash # Name of the .SWI or .PWC file
pwcName = "Koc10000"
# If there's a proper hourly precipitation file, use that, otherwise generate the precipitation event programmatically from the daily value
useHPF = false
stormLengthInHours = 8 #only read when useHPF is false

pesticideEquation = 3 # 1:Sabbagh;2:refitSabbagh;3:mass-bal.;4:Chen

# Ksat - Any positive value in m/s
# Set to any value -9999 or less to use the lookup table method (based on soil texture pulled from the PWC scenario)
# Set to any 0 > value > -9999 to use the Wösten et al., 1999 formula (based on soil texture, bulk density, and organic matter from PWC scenario, and specific to topsoils)
Ksat = -4999 #in m/s
    
# The dimensions of the VFS are calculated from its width (which VFS calls its length) and the area of the pond into which it flows
# It can be a round pond with a round VFS around it, a square pond with a square VFS around it, or a square pond with a rectangular VFS on one side
shapeFlag = 1 # 1 for round, 2 for square, 3 for rectangular


#Less often changed are the paths to unchanging files
# Executables file must contain vfsm,exe and vvwm.exe, as well as SoilTriangle.csv
exePath = "Z:\\SharedwithVM\\VFS\\executables\\" # must end in a double back-slash or slash
# The path to the standard scenarios
# If a .SCN file by the name of that used for the PWC run is found in the workingPath, that is used, otherwise a version is required in scenPath
scenarioPath = "Z:\\SharedwithVM\\Scenarios\\" # must end in a double back-slash or slash

# The location of the pre-run .zts files where θ, the pre-storm water content of the VFS, is found.
turfPath = "Z:\\SharedwithVM\\VFS\\CanadianTurfZts\\" # must end in a double back-slash or slash

# Everything below this ought not to require user input
# To allow multiple simulations in the same folder
if isfile(string(workingPath, "vvwmTransfer.txt"))
    scenarioName = readlines(string(workingPath, "vvwmTransfer.txt"))[29]
else
    println("VFSPipe has already been run on these PWC results. Please re-run PWC and try again.")
    return
end

firstRun = true
for width in stripWidthInM
    
    lastrun = false
    if width == width[end]
        lastrun = true
    end

    # All of the user input variables collected here are passed to the other module in a single struct called a userInputs
    #usInp = userInputs(projectName = projectName, stripWidthInM = width, twoCharacterCode = twoCharacterCode, workingPath = workingPath, pwcName = pwcName, useHPF = useHPF, stormLengthInHours = stormLengthInHours, exePath = exePath, turfPath = turfPath, pesticideEquation = pesticideEquation, Ksat=Ksat, shapeFlag = shapeFlag, isFirstRun = firstRun, scenarioName = scenarioName, scenarioPath = scenarioPath)
    usInp = userInputs(stripWidthInM = width, workingPath = workingPath, pwcName = pwcName, useHPF = useHPF, stormLengthInHours = stormLengthInHours, exePath = exePath, turfPath = turfPath, pesticideEquation = pesticideEquation, Ksat=Ksat, shapeFlag = shapeFlag, isFirstRun = firstRun, isLastRun = lastrun, scenarioName = scenarioName, scenarioPath = scenarioPath)

    #This runs VFSMOD for all precipitation events, rewrites the zts file and returns the string to run VVWM
    newVVWM, oldTXT, newTXT = vfsMain(usInp)
    run(newVVWM)
    noVFS = getVVWMText(columnNames(),oldTXT)
    yesVFS = getVVWMText(columnNames(),newTXT)
    plot(noVFS.Run, [noVFS.Peak,yesVFS.Peak], title = string("Efficacy of ", width, "m VFS for EEC Reduction"), label = ["Peak Without VFS" "Peak With VFS"], linecolor = ["Brown" "Green"],lw = 2, xlabel = "Simulation Year", ylabel = "PPB")
    plot!(noVFS.Run, [noVFS.OneYear,yesVFS.OneYear], title = string("Efficacy of ", width, "m VFS for EEC Reduction"), label = ["Yearly Without VFS" "Yearly With VFS"], linecolor = ["Dark Orange" "Light Green"],lw = 2, xlabel = "Simulation Year", ylabel = "PPB")
    global firstRun = false
end

#A holdover from VFSMOD which may have other applications
#User choice - Can replace the last two letters of the input files for the temporary VFSMOD working files so they do not get written over
##projectName = "VFSM" #Six characters only ### Currently not used
##twoCharacterCode = "OD" ### Currently not used