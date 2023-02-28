# This takes the .zts output from PRZM5 to parameterize the input files of VFSMOD
# for each individual runoff event, and compiles the results into
# another .zts file upon which VVWM is run
include("./VFSUtils.jl")

#User choice - an array of widths in meters
stripWidthInM = [10]

# The path to the working directory - should contain the .SWI or .PWC PWC file, and the output files from the original run (.zts, .inp, vvwmtransfer.txt).

workingPath = "C:\\Users\\hcuser\\Documents\\Local VFS\\HPF\\"

# Name of the .PWC file
pwcName = "GP-GAOnion"

# *******************************************************************************************************************************************************************
#Less often changed are the paths to unchanging files
# Executables folder must contain PRZM5.exe, vfsm.exe and vvwm.exe
exePath = "C:\\Users\\hcuser\\Documents\\Local VFS\\executables\\" # must end in a double back-slash or slash

# The location of the pre-run .zts files where θ, the water content of the VFS soil, is found.
# If blank (just "") a new .zts with the appropriate θ data will be created
# thetaPath = "Z:\\SharedwithVM\\VFS Connector Demo\\" # must end in a double back-slash or slash
thetaPath = ""

# The path to the standard scenarios is not required any longer - all scenario information is read from the .PWC or .SWI file

# *******************************************************************************************************************************************************************
# The following all have defaults, but are of interest for the model's sensitivity to them
# If true, VFSPipe will look in the weather directory for a .hpf file with the same name as the .dvf or .wea file used to run PWC
useHPF = true
stormLengthInHours = 8 # RECOMMENDED CHOICE IS 8

pesticideEquation = 3 # 1:Sabbagh;2:refitSabbagh;3:mass-bal.;4:Chen - RECOMMENDED CHOICE IS 3
remobilizationFlag = 1 # Recommended choice is 1 - only the part of the dissolved pesticide in mixing layer in equilibrium with moisture content is remobilized at the start of the next event

# Ksat - Any positive value in m/s
# Set to any value -9999 or less to use the lookup table method (based on soil texture pulled from the PWC scenario)
# Set to any 0 > value > -9999 to use the Wösten et al., 1999 formula (based on soil texture, bulk density, and organic matter from PWC scenario, and specific to topsoils)
Ksat = -9999 #in m/s - RECOMMENDED CHOICE IS -9999

# The dimensions of the VFS are calculated from its width (which VFSMOD calls its length) and the area of the pond into which it flows
# The default is a square field based on the field area with a rectangular VFS on one edge
# It can also be a round pond with a round VFS around it, a square pond with a square VFS around it,
# ...a rectangular pond with a rectangular VFS on one side (based on Hydraulic Length from PWC), or a square pond based on the pond area with a rectangular VFS and field on one side
shapeFlag = 1 # 1 for default, 2 for round, 3 for square, 4 for rectangular with a length determined from the hydraulic length and area of the field, anything else for rectangular on one side of a square pond

# Slope - an array of any positive or zero values in percent 
# Any negative value to read slope from PWC scenario


slopesInPercent = [-9999]



# Everything below this ought not to require user input

#***********************************************************************************************************************************************************************
# A check for a 'fresh' vvwmTransfer file - if it exists, capture the scenario name to use as a unique name for the working files
# Unique names allow multiple, consecutive simulations in the same folder
oldTXT = ""
newTXT = ""
if !isfile(string(workingPath, "vvwmTransfer.txt"))
    #scenarioName = readlines(string(workingPath, "vvwmTransfer.txt"))[29]
    error("VFSPipe has already been run on these PWC results. Please re-run PWC and try again.")
end

firstRun = true
for width in stripWidthInM

    for slope in slopesInPercent

        lastrun = false
        if width == width[end]
            lastrun = true
        end

        # All of the user input variables collected here are passed to the other module in a single struct called a userInputs
        usInp = userInputs(stripWidthInM=width, workingPath=workingPath, pwcName=pwcName, useHPF=useHPF, stormLengthInHours=stormLengthInHours, exePath=exePath, thetaPath=thetaPath, pesticideEquation=pesticideEquation, Ksat=Ksat, shapeFlag=shapeFlag, slope=slope, isFirstRun=firstRun, isLastRun=lastrun, remobilizationFlag=remobilizationFlag)#, scenarioName = scenarioName, scenarioPath = scenarioPath)

        #This runs VFSMOD for all runoff events, rewrites the zts file and returns the string to run VVWM
        newVVWM, oTXT, nTXT = vfsMain(usInp)
        run(newVVWM)

        # multiple assignments don't seem to work with global variables, so we do it locally, then globally
        global oldTXT = oTXT
        global newTXT = nTXT
        global firstRun = false
    end
end

#A holdover from VFSMOD which may have other applications
#User choice - Can replace the last two letters of the input files for the temporary VFSMOD working files so they do not get written over
##projectName = "VFSM" #Six characters only ### Currently not used
##twoCharacterCode = "OD" ### Currently not used

# This will only work on the last run, but was requested by users
# This is obviously going to be improved upon, adding the ability to 'browse' between results
 noVFS = getVVWMText(columnNames(),oldTXT)
 yesVFS = getVVWMText(columnNames(),newTXT)
 plot(noVFS.Run, [noVFS.Peak,yesVFS.Peak], title = string("Efficacy of ", stripWidthInM[end], "m VFS for EEC Reduction"), label = ["Peak Without VFS" "Peak With VFS"], linecolor = ["Brown" "Green"],lw = 2, xlabel = "Simulation Year", ylabel = "PPB")
 plot!(noVFS.Run, [noVFS.OneYear,yesVFS.OneYear], title = string("Efficacy of ", stripWidthInM[end], "m VFS for EEC Reduction"), label = ["Yearly Without VFS" "Yearly With VFS"], linecolor = ["Dark Orange" "Light Green"],lw = 2, xlabel = "Simulation Year", ylabel = "PPB")
