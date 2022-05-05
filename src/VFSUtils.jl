### This file contains all of the structs and the functions that manipulate them
# that are needed to run VFSMOD from PRZM outputs, and VVWM from VFSMOD outputs

# Notes marked "2022 Guidance" refer to Ritter and Muñoz-Carpena, VFSMOD Input Parameter Guidance for Long-tern Regulatory Scenarios

using CSVFiles, DataFrames, PrettyTables, Printf, CSV, Parameters, Plots, StaticArrays #Gtk removed

# THIS SHOULD BE REPLACED BY IMPLEMENTING THE "UNITS" PACKAGE
#constants, just to avoid magic numbers
@with_kw struct constants
    @deftype Int64
    cmInAMeter = 100
    secondsInAnHour = 3600
    hoursInADay = 24
    mSqInAHa = 10000
    gInATonne = 1000000
    cm3InAM3 = 1000000
    mgInAGram = 1000
    expectedRefTempInC = 20 #VFSMOD expects the degradation rate at this temperature
    ρQuartz::Float64 = 2.65 #Density of quartz in tonnes per m^3 or g cm^3
    ρWater::Float64 = 1.00 #Density of Water in g per cm^3
    cent = 100 #As in percent
    kgInATonne = 1000
    cmSqInAMSq = 10000
end

# matches VFSMOD parameter names
@with_kw struct soilParameters{R<:Real}
    #Using the variable names from VFSMOD
    θI::R = 0.0 #This is the only value that changes from run to run
    SM::R = 0.0 # Maximum surface storage
    SCHK::R = 0.5 # relative distance from upper edge..?
end

# matches VFSMOD parameter names
@with_kw struct chemicalParameters{R<:Real}
    IKD::Int8 = 1 #1=KOC used
    KOCD::R = 100.0 #
    DGGHALF::R = 30.0 #Halflife, not K
    DGML::R = 0.05 #dispersion length of chemical, meters, default from 2022 guidance
end

#These are the parameters required by the water quality file, and are NOT sourced from other structs
# Other structs required for the output include chemical and scenario parameters
@with_kw struct waterQualityParameters{R<:Real}
    IQPRO::Int8 = 3 # 1:Sabbagh;2:refitSabbagh;3:mass-bal.;4:Chen
    IDG::Int8 = 2 #2=degradation by decay rate only
    dgML::R = 2.0 #Surface mixing layer thickness = constant
    IMOB::R = 1.0 # Portion of the residue available for carryover the the next runoff event
end

# matches VFSMOD parameter names
@with_kw struct sedimentParameters{R<:Real}
    #Using the variable names from VFSMOD
    NPART::Int8 = 7 # 7 = read user input values
    CI::R = 0.0 #This is the only value that changes from run to run
    SG::R = 2.65 #This is the density of quartz
end

# A blend of PWC and VFSMOD values
@with_kw struct scenarioParameters{R<:Real}
    # All of these 'defaults' are overwritten by values from the scenario (.SCN), or calculated from information in the scenario

    #For a few calculations
    fieldAreaInHa::R = 10.0 #in ha
    pondAreaInM2::R = 10000.0 #in m^2
    hydraulicLengthInM::R = 356.8 #in m
    #For water quality
    OCP::R = 3.2 #Organic Carbon %
    CCP::R = 19.0 #Clay %
    FC::R = 0.31 # Field capacity
    # For sediments
    COARSE::R = 0.75 # really just sand fraction
    POR::R = 0.472 # same as saturated water content
    DP::R = 0.0098 #This should be modified for each SCENARIO
    #For Soil
    VKS::R = 6.06E-6 # Saturated hydraulic conductivity
    SAV::R = 0.1101 # Green-Ampts' suction at wet front
    #Note in appendix to VFSMOD Manual says 'equivent to POR value'
    θSoil::R = 0.47 # Saturated soil-water content (same as porosity for now)
    ρSed::R = 2.65 # Density of sediment particles
    # Slope needed for filterParameters
    SOA::R = 0.0200 # As a fraction
end

# matches PRZM parameter names
@with_kw struct waterQualityOutput
    @deftype Float64
    RUNF0 = 0.0 #cm per ha per day - water per field area
    ESLS0 = 0.0 #tonnes per day - solids per field area
    RFLX1 = 0.0 # g per sqcm per day - mass of pesticide by runoff
    EFLX1 = 0.0 # g per sqcm per day - mass of pesticide by erosion
    DCON1 = 0.0 # g per sqcm per day - mass of degradate 1 by runoff
    INFL0 = 0.0 # g per sqcm per day - mass of degradate 1 by erosion
end
# matches VFSMOD parameter names
@with_kw struct filterParameters
    @deftype Float64
    FWIDTH = 363.821 # Width of strip in m
    VL = 3.0 # Length of strip in m
    N::Int64 = 11  # Number of nodes
    θw = 0.5 # Time weighted factor, Crank-Nicholson solution
    CR = 0.8 # Courant number
    MAXITER::Int64 = 350 # Maximum number of interations
    NPOL::Int64 = 3 # Number of nodal points
    IELOUT::Int64 = 1 #Flag for output
    KPG::Int64 = 1 # Flag for Petrox-Galerkin solution
    NPROP::Int64 = 1 #Number of segments
    SX = 3.0 #Distance of filter of uniform surface properties - default is same as VL
    RNA = 0.45 #Manning's n for segment
    SOA = 0.0200 # Slope not read here - write method takes from scenarioParameters
    IWQ::Int64 = 1 #Flag for water quality calculation, 1 = do it
    areaInMSq = 1091.5 # Area in m^2
    areaInHa = 0.1091 # Area in Ha
end

# matches VFSMOD parameter names
@with_kw struct grassParameters
    @deftype Float64
    SS = 2.15
    VN = 0.012
    H = 18
    VN2 = 0.05
    ICO::Int8 = 0
end

@with_kw struct rainInput
    nrain::Int64 = 0
    rpeak::Float64 = 0.0
    rain::SArray{Tuple{200,2},Float64}
end

@with_kw struct runoffInput
    @deftype Float64
    swidth = 0.0
    slength = 0.0
    nbcroff::Int64 = 0
    bcropeak = 0.0
    bcroff::SArray{Tuple{200,2},Float64}
end

# A container for all of the filenames and some paths required to make this crazy thing go
# Easy and safe and portable
@with_kw struct inOutFileNames
    @deftype String

    przmInName = ""
    vfsmod::Cmd = ""
    vvwm::Cmd = ""
    thetaInName = ""
    przmOutName = ""
    hpfInName = ""
    prjInName = ""
    prjOutName = ""
    scnFileName = ""
    swiFileName = ""
    isPWCFormat::Bool = false
    vvwmTransferFileName = ""
    dailyWeatherFileName = ""
    ikwOutName = ""
    irnOutName = ""
    igrOutName = ""
    iroOutName = ""
    isoOutName = ""
    isdOutName = ""
    iwqOutName = ""
    owqFileName = ""
    osmFileName = ""
end

# A container for all of the column headers
@with_kw struct columnNames
    @deftype Array{String}
    ztsColumns = ["Yr", "Mo", "Dy", "RUNF0", "ESLS0", "RFLX1", "EFLX1", "DCON1", "INFL0"]
    hpfColumns = ["Date", "Duration", "Total", "Midnight", "OneAM", "TwoAM", "ThreeAM", "FourAM", "FiveAM", "SixAM", "SevenAM", "EightAM", "NineAM", "TenAM", "ElevenAM", "Noon", "OnePM", "TwoPM", "ThreePM", "FourPM", "FivePM", "SixPM", "SevenPM", "EightPM", "NinePM", "TenPM", "ElevenPM"]
    dvfColumns = ["Date", "Total", "A", "Temperature", "B", "C"]
    dvf14Columns = ["Date", "Total", "A", "Temperature", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"]
    weaColumns = ["Month", "Day", "Year", "Total", "A,", "Temperature", "B", "C"]
    vvwmTxtColumns = ["Run", "Peak", "FourDay", "TwentyOneDay", "SixtyDay", "NinetyDay", "OneYear", "PWPeak", "PWTwentyOneDay"]
    vvwmCsvColumns = ["Depth", "Average", "PWAverage", "Peak"]
end
# A container for all the things that the user can change easily
@with_kw struct userInputs
    @deftype String
    stripWidthInM::Float64 = 0.0
    workingPath = ""
    pwcName = ""
    useHPF::Bool = false
    stormLengthInHours::Int64 = 8
    exePath = ""
    thetaPath = ""
    pesticideEquation::Int64 = 3 # 1:Sabbagh;2:refitSabbagh;3:mass-bal.;4:Chen
    Ksat::Float64 = 0.0
    shapeFlag::Int64 = 1
    slope::Float64 = 0.0 # In percent
    isFirstRun::Bool = true
    isLastRun::Bool = false
    scenarioName = ""
    scenarioPath = ""
end

# Generate the table of soil classes here, so it doesn't have to be read in from some external file
function soilClasses()
    scs = [[0, 45, 0, 40, 40, 100, "Clay", 0.000000167, 3.16E-01, 0.0002], [0, 20, 40, 60, 40, 60, "SiltyClay", 0.000000278, 2.92E-01, 0.0024], [45, 65, 0, 20, 35, 55, "SandyClay", 0.000000333, 2.39E-01, 0.0066], [0, 20, 40, 73, 27, 40, "SiltyClayLoam", 0.000000556, 2.73E-01, 0.0025], [45, 80, 0, 28, 20, 35, "SandyClayLoam", 0.000000833, 2.19E-01, 0.0091], [20, 45, 15, 52, 27, 40, "ClayLoam", 0.000000556, 2.09E-01, 0.0018], [0, 20, 88, 100, 0, 12, "Silt", 0, 0.00E+00, 0.0019], [0, 50, 50, 88, 0, 27, "SiltyLoam", 0.00000189, 1.67E-01, 0.0027], [23, 52, 28, 50, 7, 27, "Loam", 0.00000367, 8.89E-02, 0.0035], [50, 70, 0, 50, 0, 20, "SandyLoam", 0.00000606, 1.10E-01, 0.0098], [70, 86, 0, 30, 0, 15, "LoamySand", 0.0000166, 6.13E-02, 0.0135], [86, 100, 0, 14, 0, 10, "Sand", 0.0000654, 4.95E-02, 0.017]]
    scs = DataFrame((SandLow=r[1], SandHigh=r[2], SiltLow=r[3], SiltHigh=r[4], ClayLow=r[5], ClayHigh=r[6], Class=r[7], Ksat=r[8], Sav=r[9], DP=r[10]) for r in scs)
    return scs
end
#The .zts file uses an inconvenient format for the numbers, required for the correct spacing
#This function just formats numbers that way - in base 10 scientific notation with one digit left of the
#decimal, 4 digits right of the decimal, and three digits in the exponent
#It turns out that this doesn't matter, but it looks 'right'
function likePrzm(num)
    return string(@sprintf("%.4E", num)[1:8], @sprintf("%03d", parse(Int64, @sprintf("%.4E", num)[9:10])))
end

#This function generates all of the filenames used by VVWM and VFSMOD
# It takes most of the user input...
#It writes the project (.prj) file required by VFSMOD
# IFF this is the first run of a batch, it makes a copy of the vvwmtransfer.txt file, to avoid mistakes while running several instances in the same folder
#It returns an inOutNames struct

function writeFileNames(ui::userInputs)
    # if length(ui.projectName) > 6
    #     error("Project name must be <= 6 characters")
    # end
    projectName = "VFSM"

    # Append a run number to allow running more than one instance in the same directory simultaneously
    runNumber = 0
    while isfile(string(ui.workingPath, projectName, runNumber, ".igr"))
        runNumber += 1
    end

    executeString = string(projectName, runNumber, ".prj")
    exePath = ui.exePath[1:end-1] #To make the path 'command friedly'
    vfsmod = `$exePath\\vfsm.exe $executeString 0` #The zero prevents the textual output from VFSMOD to the REPL

    # set the PWC file name from the user input
    # if the user has not specified the file extension, try appending .SWI or .PWC
    if isfile(string(ui.workingPath, ui.pwcName))
        swiFileName = string(ui.workingPath, ui.pwcName)
    else
        swiFileName = string(ui.workingPath, ui.pwcName, ".SWI")
        if !isfile(swiFileName)
            swiFileName = string(ui.workingPath, ui.pwcName, ".PWC")
        end
        if !isfile(swiFileName)
            error("Unable to find a .SWI or .PWC file with the name ", ui.pwcName)
        end
    end

    # The PWC2 format is slightly different than the PWC1.52 format
    # Check the end of the first line of the file for the model version number
    if isfile(swiFileName)
        swi = open(swiFileName, "r")
    else
        error(string("Unable to find ", swiFileName))
    end
    isPWCFormat = parse(Float64, split(readline(swi))[end]) >= 2
    close(swi)

    vvwmTransferFileName = string(ui.workingPath, ui.pwcName, "_", ui.scenarioName, "_vvwmTransfer.txt")
    # IF this is the first run after a PRZM run, make a safe copy - allows direct comparison of what enters the pond
    if ui.isFirstRun
        cp(string(ui.workingPath, "vvwmTransfer.txt"), string(vvwmTransferFileName[1:end-4], "_sansVFS.txt"), force=true)
        rm(string(ui.workingPath, "vvwmTransfer.txt"))
    end

    vvwm = `$exePath\\VVWM.exe $vvwmTransferFileName`

    # For debugging convenience
    #if !isPWC2Format
    #    vvwm = `$exePath\\VVWM-forPWC1.52.exe $vvwmTransferFileName`
    #end

    przmInName = string(ui.workingPath, projectName, "_before.zts")
    #przmInName = string(workingPath, projectName, "NM.zts")
    cp(string(ui.workingPath, ui.pwcName, ".zts"), przmInName, force=true)
    przmOutName = string(ui.workingPath, ui.stripWidthInM, "m_", projectName, ".zts")
    #For now, since we don't really need an HPF - could create from dvf or met right here
    #In future, should explore real hourly with the 'new' met data
    hpfInName = string(ui.workingPath, projectName, ".HPF")
    prjInName = string(ui.workingPath, projectName, ".prj")
    prjOutName = string(ui.workingPath, projectName, runNumber, ".prj")

    #scnFileName = string(ui.workingPath, readlines(swiFileName)[49], ".SCN")
    # No longer need anything from the .SCN or .SCN2 file
    #scnFileName = string(ui.workingPath, ui.scenarioName, ".SCN") #Use the scenario file name from the VVWMTransfer file.
    #if isfile(scnFileName)
    #    # Do nothing
    #else
    #    scnFileName = string(ui.workingPath, ui.scenarioName, ".SCN2")
    #    if isfile(scnFileName)
    #        # Do nothing
    #    else        
    #        scnFileName = string(ui.scenarioPath, readlines(swiFileName)[49], ".SCN")
    #        if isfile(scnFileName)
    #            #Do nothing
    #        else
    #            scnFileName = string(ui.scenarioPath, readlines(swiFileName)[49], ".SCN2")
    #        end
    #    end
    #end

    dailyWeatherFileName = readlines(swiFileName)[50]
    if isfile(dailyWeatherFileName)
        # Do nothing
    else
        dailyWeatherFileName = string(ui.workingPath, split(readlines(scnFileName)[2], "\\")[end])
    end
    thetaInName = string(ui.thetaPath, split(readlines(swiFileName)[50], "\\")[end][1:(end-4)], "Theta.zts")

    #Create a project file for VFSMOD to read the file names too
    fileTypes = ["ikw", "iso", "igr", "isd", "irn", "iro", "iwq", "og1", "og2", "ohy", "osm", "osp", "owq"]
    prj = open(prjOutName, "w")
    for i in 1:length(fileTypes)
        write(prj, string(fileTypes[i], "=", projectName, runNumber, ".", fileTypes[i]), "\n")
    end
    close(prj)

    # The input files for VFSMOD that are changed by this program
    ikwOutName = readlines(prjOutName)[1][5:end]
    irnOutName = readlines(prjOutName)[5][5:end]
    iroOutName = readlines(prjOutName)[6][5:end]
    isoOutName = readlines(prjOutName)[2][5:end]
    isdOutName = readlines(prjOutName)[4][5:end]
    iwqOutName = readlines(prjOutName)[7][5:end]
    igrOutName = readlines(prjOutName)[3][5:end]

    # The results of the VFSMOD run have two purposes:
    # To create the output for VVWM
    # And to create the input for the next run of VFSMOD
    owqFileName = string(ui.workingPath, readlines(prjOutName)[13][5:end])
    osmFileName = string(ui.workingPath, readlines(prjOutName)[11][5:end])

    #Place all of these in a struct so they cannot be accidentally overwritten
    # scnFileName is removed since all pertinent information is in the .SWI or .PWC file
    inOutNames = inOutFileNames(vfsmod=vfsmod, vvwm=vvwm, przmInName=przmInName, thetaInName=thetaInName, przmOutName=przmOutName, hpfInName=hpfInName, prjInName=prjInName, prjOutName=prjOutName, swiFileName=swiFileName, isPWCFormat=isPWCFormat, dailyWeatherFileName=dailyWeatherFileName, ikwOutName=ikwOutName, irnOutName=irnOutName, iroOutName=iroOutName, isoOutName=isoOutName, isdOutName=isdOutName, iwqOutName=iwqOutName, owqFileName=owqFileName, igrOutName=igrOutName, osmFileName=osmFileName, vvwmTransferFileName=vvwmTransferFileName)
    return inOutNames
end

function checkFilesExist(fns::inOutFileNames, ui::userInputs)
    missingFiles = []
    if !ui.useHPF
        inNames = [fns.przmInName, fns.thetaInName, fns.swiFileName, fns.dailyWeatherFileName]
    else
        inNames = [fns.przmInName, fns.thetaInName, fns.hpfInName, fns.swiFileName, fns.dailyWeatherFileName]
    end
    for fn in inNames
        if !isfile(fn)
            push!(missingFiles, fn)
        end
    end
    return missingFiles
end

# This function writes the .irn file required by VFSMOD
# There are two versions: one which interprets daily precipitation files (for all Canadian scenarios) and one which interprets hourly precipitation files
#This takes a vector of HOURLY rain intensities and writes out an irn file of durations and intensities
#This version trims zeroes off both ends to reduce run time of VFSMOD
function writePrecipitation(row::DataFrame, fns::inOutFileNames)
    #Containers
    times = Vector{Int16}()
    rates = Vector{Float32}()
    #Convert this row from a DataFrame into a vector for easy manipulation
    row = vec(permutedims(Vector(row)))
    #Remove any zeroes from the beginning
    while row[1] == 0
        deleteat!(row, 1)
    end
    #Remove any zeroes from the end
    while row[end] == 0
        deleteat!(row, length(row))
    end

    #Do not want reduduncancy, but easier to write them all and remove than to control what to write
    for i in 1:length(row)
        t = (i - 1) * constants().secondsInAnHour # We want the hour, not the index
        times = [times; t]
        r = row[i] / constants().cmInAMeter / constants().secondsInAnHour #HPF is in cm/h, irn is in m/s
        rates = [rates; r]
    end
    i = 1 #length(times)
    while length(times[i:end]) > 2
        if rates[i] == rates[i+1] #&& rates[i] == rates[i+2]
            deleteat!(times, i + 1)
            deleteat!(rates, i + 1)
        else
            i += 1
        end
    end

    #The runoff calculation needs the precipitation time
    #Need to record it before adding on the zeroes
    stormTime = maximum(times)

    #According to the UFSMOD manual, it is generally it is not raining at the end of the simulation
    #The example pipe appears to end the simulation with a two-hour rain-less period
    if maximum(times) < (constants().hoursInADay - 2) * constants().secondsInAnHour
        times = [times; (maximum(times) + 1)] #one second after the last non-zero
        rates = [rates; 0]
        times = [times; ((maximum(times) - 1) + 2 * constants().secondsInAnHour)] #two hours after the last non-zero
        rates = [rates; 0]
    end

    # create the struct before altering the vectors for text output
    rain = zeros(Float64, 200, 2)
    rain[1:length(times), 1:2] = hcat([times], [rates])
    storm = rainInput(nrain=numberOfDataInAFakeStorm, rpeak=secondlyRate, rain=rain)

    #VFSMOD expects a description of the data at the top of the data
    pushfirst!(times, length(times))
    pushfirst!(rates, maximum(rates))
    irn = open(fns.irnOutName, "w")
    for j in 1:length(times)
        write(irn, string(@sprintf("%9s", times[j]), "  ", @sprintf("%.4E", rates[j])), "\n")
        #write(irn, string(" ", @sprintf("%08u",times[j]), "  ", @sprintf("%.4E",rates[j])), "\n")
    end
    close(irn)
    return stormTime, storm
end
# This version takes a daily rain amount and generates a stormTime seconds long storm from it
# It has been updated to reflect the 2022 guidance from the VFS Working Group
# The hyetograph is rectangular
function writePrecipitation(precip::Float64, stormTime::Int64, fns::inOutFileNames)

    numberOfDataInAFakeStorm = 4
    startOfStorm = 0
    oneSecondAfterStorm = stormTime + 1
    secondlyRate = precip / constants().cmInAMeter / stormTime #.dvf is in cm/storm, rainfall intensity in .irn is in m/s
    # There are two hours of simulation after the storm
    # endOfSimulationTime = stormTime + 2*constants().secondsInAnHour
    # 2022 guidance suggests 1.25 times the storm length instead of 2 hours (which is 2 hours for an 8-hour storm)
    endOfSimulationTime = stormTime * 1.25
    # Unless, of course, the storm lasts for > 22 hours
    if endOfSimulationTime > constants().hoursInADay * constants().secondsInAnHour
        endOfSimulationTime = constants().hoursInADay * constants().secondsInAnHour
    end

    irn = open(fns.irnOutName, "w")
    # The first line has the number of events to read and the maximum rate, which is the same as the rate, since the hydrograph is square
    write(irn, string(@sprintf("%9s", numberOfDataInAFakeStorm), "  ", @sprintf("%.4E", secondlyRate)), "\n")
    write(irn, string(@sprintf("%9s", startOfStorm), "  ", @sprintf("%.4E", secondlyRate)), "\n")
    write(irn, string(@sprintf("%9s", stormTime), "  ", @sprintf("%.4E", secondlyRate)), "\n")
    write(irn, string(@sprintf("%9s", oneSecondAfterStorm), "  ", @sprintf("%.4E", 0.0)), "\n")
    write(irn, string(@sprintf("%9s", endOfSimulationTime), "  ", @sprintf("%.4E", 0.0)), "\n")
    close(irn)
    #rain = Array{Float64,2}(undef,200,2)
    rain = zeros(Float64, 200, 2)
    rain[1:numberOfDataInAFakeStorm, 1:2] = hcat([startOfStorm, stormTime, oneSecondAfterStorm, endOfSimulationTime], [secondlyRate, secondlyRate, 0.0, 0.0])
    storm = rainInput(nrain=numberOfDataInAFakeStorm, rpeak=secondlyRate, rain=rain)
    return storm
end

# this takes the scenarioParameters struct and writes out the .iso file for VFSMOD
function writeSoil(θ, scen::scenarioParameters, fns::inOutFileNames)
    soil = soilParameters(θI=θ)
    iso = open(fns.isoOutName, "w")
    write(iso, string("  ", @sprintf("%.2E", scen.VKS), "  ", @sprintf("%.4f", scen.SAV), "  ", @sprintf("%.4f", scen.θSoil), "  ", @sprintf("%.4f", soil.θI), "  ", @sprintf("%.4f", soil.SM), "  ", @sprintf("%.4f", soil.SCHK)))
    close(iso)
end

#This function takes the daily runoff from PRZM, a storm length, and a filemane
#and writes an iro file
#This has been updated to reflect the 2022 guidance from the VFS Working Group
function writeRunoff(RUNF0, scen::scenarioParameters, filter::filterParameters, eventLength, fns::inOutFileNames)

    fieldAreaInMSq = scen.fieldAreaInHa * constants().mSqInAHa #in ha - confirmed in VVWM manual that runoff is normalized to FIELD area (not waterbody area)

    SWIDTH = filter.FWIDTH#(sqrt(filter.areaInMSq/π)+filter.FWIDTH)*2π
    SLENGTH = fieldAreaInMSq / SWIDTH

    #timeStep::Int64 = (constants().secondsInAnHour/2) # Arbitrary choice - half an hour
    #The 8 hour storms in the example have 10 time steps - it is not known why
    #For fun, here, the time step is set and the number derivied - which will allow the timestep to be a user input in future
    #numberOfSteps::Int16 = eventLength/timeStep
    NBCROFF::Int64 = 3 # As per 2022 guidance for a triangular hydrograph
    #Get from cm/ha/day to m^3/day, and then per second
    runOn = RUNF0 / constants().cmInAMeter * fieldAreaInMSq / eventLength
    BCROPEAK = 2 * runOn
    #Containers
    times = Vector{Int16}()
    rates = Vector{Float32}()

    #write the times and the rate for the event
    #=     for i in 0:numberOfSteps
            push!(times,i*timeStep)
            push!(rates,runOn)
        end =#

    # Write the times and rate for the triangular hydrograph
    push!(times, 0)
    push!(rates, 0)
    push!(times, Int(floor(eventLength / 2.67 + 0.5)))
    push!(rates, BCROPEAK)
    push!(times, eventLength)
    push!(rates, 0)

    #The second line has the number of timesteps and the peak rate of the hydrograph
    pushfirst!(times, NBCROFF)
    pushfirst!(rates, BCROPEAK)

    iro = open(fns.iroOutName, "w")
    write(iro, string(" ", SWIDTH, " ", SLENGTH), "\n")
    for i = 1:length(times)
        write(iro, string(@sprintf("%9s", times[i]), " ", @sprintf("%.4E", rates[i])), "\n")
    end
    close(iro)
end

# this populates the scenarioParameters struct from the information in the PWC .SCN File
# It causes more problems than it's worth another version reads from the .SWI or .PWC directly
function readScenarioParameters(scnFileName, usInp::userInputs)
    scn = open(scnFileName, "r")
    #Throw away the first three lines
    for i in 1:6
        readline(scn)
    end
    fieldArea = 10.0
    fldArea = readline(scn) #7
    if fldArea != ""
        fieldArea = parse(Float64, fldArea) / constants().mSqInAHa
    end
    pondArea = 10000.0
    pndArea = readline(scn) #8
    if pndArea != ""
        pondArea = parse(Float64, pndArea)
    end
    for i in 9:13
        readline(scn)
    end
    benthicPorosity = readfirst(scn) # 14
    ρBenthicBulk = readfirst(scn) # 15
    for i in 16:48
        readline(scn)
    end
    slope = readfirst(scn) / constants().cent #49 - must be fractional, not percent
    slpPrct::Float64 = usInp.slope
    if slpPrct >= 0 # The user has specified a slope (in percent) to override the scenario
        slope = slpPrct / constants().cent # Must be a fraction, not a percent
    end
    hydraulicLength = readfirst(scn) #50 in m
    for i in 51:52
        readline(scn)
    end
    ρ = readfirst(scn) #53 in g per cm^3 (or tonnes per m^3)
    fieldCapacity = readfirst(scn) #54
    wiltingPoint = readfirst(scn) #55
    ocPercent = readfirst(scn) #56
    readline(scn) #57
    sandPercent = readfirst(scn) #58
    clayPercent = readfirst(scn) #59
    close(scn)
    POR = 1 - (ρ / constants().ρQuartz)
    texture = getSoilClass(sandPercent, (100 - sum([sandPercent, clayPercent])), clayPercent)
    #This derivation probably highlights a logical flaw in determining the bulk density of sediments
    #In PWC scenario development rather than providing a better value for sediment bulk density
    #It shall be left here but not used in the writing of the sediment input file for now
    ρSedimentParticle = (ρBenthicBulk - benthicPorosity * constants().ρWater) / (1 - benthicPorosity)

    #Alternate methods to determine Ksat (a.k.a VKS)
    topsoil = 1 #This is a qualitative binary in the equation - set here to 1 so that all VFS soils are topsoils when using Wösten
    VKS::Float64 = usInp.Ksat
    if VKS < -9998.999 #it's not an integer, so just being cautious
        VKS = texture[1]
    elseif -9999 < VKS < 0.0 # use the empirical equation from Wösten et al., 1999 (Geoderma 90 1999 169–185)
        silt = (100 - sum([sandPercent, clayPercent]))
        om = ocPercent * 1.724 # Empirical adjustment from organic carbon percent to organic matter percent
        Wösten = exp(7.755
                     + 0.0352 * silt
                     + 0.93 * topsoil
                     -
                     0.967 * ρ^2
                     -
                     0.000484 * clayPercent^2
                     -
                     0.000322 * silt^2
                     +
                     0.001 / silt
                     -
                     0.0748 / om
                     -
                     0.643 * log(silt)
                     -
                     0.01398 * ρ * clayPercent
                     -
                     0.1673 * ρ * om
                     +
                     0.02986 * topsoil * clayPercent
                     -
                     0.03305 * topsoil * silt)
        VKS = Wösten / constants().cmInAMeter / constants().hoursInADay / constants().secondsInAnHour
    else
        #Do Nothing: keep the user's input Ksat
    end
    scenStruct = scenarioParameters(fieldAreaInHa=fieldArea, pondAreaInM2=pondArea, hydraulicLengthInM=hydraulicLength, OCP=ocPercent, CCP=clayPercent, FC=fieldCapacity, COARSE=sandPercent / constants().cent, POR=POR, VKS=VKS, SAV=texture[2], DP=texture[3], θSoil=POR, SOA=slope)

    return scenStruct
end
# this populates the scenarioParameters struct from the information in the .SWI or .PWC file
function readScenarioParameters(fns::inOutFileNames, usInp::userInputs)
    swi = open(fns.swiFileName, "r")

    # The PWC format is slightly different than the SWI format

    # The first 75 lines are the same for both PWC and SWI
    #Throw away the first 54 lines
    for i in 1:54
        readline(swi)
    end
    fieldArea = 10.0
    fldArea = readline(swi) #55
    if fldArea != ""
        fieldArea = parse(Float64, fldArea) / constants().mSqInAHa
    end
    pondArea = 10000.0
    pndArea = readline(swi) #56
    if pndArea != ""
        pondArea = parse(Float64, pndArea)
    end
    for i in 57:61
        readline(swi)
    end
    benthicPorosity = readfirst(swi) # 62
    ρBenthicBulk = readfirst(swi) # 63

    # The rest of the values are from different locations in the SWI or PWC file
    if !fns.isPWCFormat
        for i in 64:96
            readline(swi)
        end
        slope = readfirst(swi) #97 currently in percent

        hydraulicLength = readfirst(swi) #98 in m
        for i in 99:100
            readline(swi)
        end
        ρ = readfirst(swi) #101 in g per cm^3 (or tonnes per m^3)
        fieldCapacity = readfirst(swi) #102
        wiltingPoint = readfirst(swi) #103
        ocPercent = readfirst(swi) #104
        readline(swi) #105
        sandPercent = readfirst(swi) #106
        clayPercent = readfirst(swi) #107
    else
        pwcLine = ""
        # Skip to the Soil Information section
        while pwcLine != "*** Soil Information ***"
            pwcLine = readline(swi)
        end
        readline(swi) # skip a line
        #thisLine = readline(swi)
        slope, hydraulicLength = split(readline(swi), ",")[2:3]
        slope = parse(Float64, slope)
        hydraulicLength = parse(Float64, hydraulicLength)
        while pwcLine != "*** Horizon Info *******"
            pwcLine = readline(swi)
        end
        readline(swi) # skip a line
        readline(swi) # skip a line
        ρ = readfirst(swi)
        fieldCapacity = readfirst(swi)
        wiltingPoint = readfirst(swi)
        ocPercent = readfirst(swi)
        readline(swi) # skip a line
        sandPercent = readfirst(swi)
        clayPercent = readfirst(swi)
    end

    close(swi)

    userSlope::Float64 = usInp.slope
    if userSlope >= 0 # The user has specified a slope (in percent) to override the scenario
        slope = userSlope
    end

    slope = slope / constants().cent # Must be a fraction, not a percent

    POR = 1 - (ρ / constants().ρQuartz)
    texture = getSoilClass(sandPercent, (100 - sum([sandPercent, clayPercent])), clayPercent)
    #This derivation probably highlights a logical flaw in determining the bulk density of sediments
    #In PWC scenario development rather than providing a better value for sediment bulk density
    #It shall be left here but not used in the writing of the sediment input file for now
    ρSedimentParticle = (ρBenthicBulk - benthicPorosity * constants().ρWater) / (1 - benthicPorosity)

    #Alternate methods to determine Ksat (a.k.a VKS)
    topsoil = 1 #This is a qualitative binary in the equation - set here to 1 so that all VFS soils are topsoils when using Wösten
    VKS::Float64 = usInp.Ksat
    if VKS < -9998.999 #it's not an integer, so just being cautious
        VKS = texture[1]
    elseif -9999 < VKS < 0.0 # use the empirical equation from Wösten et al., 1999 (Geoderma 90 1999 169–185)
        silt = (100 - sum([sandPercent, clayPercent]))
        om = ocPercent * 1.724 # Empirical adjustment from organic carbon percent to organic matter percent
        Wösten = exp(7.755
                     + 0.0352 * silt
                     + 0.93 * topsoil
                     -
                     0.967 * ρ^2
                     -
                     0.000484 * clayPercent^2
                     -
                     0.000322 * silt^2
                     +
                     0.001 / silt
                     -
                     0.0748 / om
                     -
                     0.643 * log(silt)
                     -
                     0.01398 * ρ * clayPercent
                     -
                     0.1673 * ρ * om
                     +
                     0.02986 * topsoil * clayPercent
                     -
                     0.03305 * topsoil * silt)
        VKS = Wösten / constants().cmInAMeter / constants().hoursInADay / constants().secondsInAnHour
    else
        #Do Nothing: keep the user's input Ksat
    end
    scenStruct = scenarioParameters(fieldAreaInHa=fieldArea, pondAreaInM2=pondArea, hydraulicLengthInM=hydraulicLength, OCP=ocPercent, CCP=clayPercent, FC=fieldCapacity, COARSE=sandPercent / constants().cent, POR=POR, VKS=VKS, SAV=texture[2], DP=texture[3], θSoil=POR, SOA=slope)

    return scenStruct
end

# Determines if b is between a and c
#Type not declared, so it works on strings too
function isBetween(b, a, c)
    if a <= b <= c
        return true
    else
        return false
    end
end

#This function takes percentages of sand, silt, and clay in that order, determines the class
#And returns a vector of Ksat, Sav, and DP in that order
function getSoilClass(sand::Real, silt::Real, clay::Real, useInternal=true)
    # Data from http://www.fao.org/tempref/FI/CDrom/FAO_Training/FAO_Training/General/x6706e/x6706e06.htm
    # Actually, that data has an error - for Silty Loam the range of Silt was wrong
    # table 4
    #Check that the values add to 100, but with some leeway for rounding error
    if isBetween(sum([sand, silt, clay]), 99.5, 100.5)
        if useInternal
            triangle = soilClasses()
        else
            triangle = DataFrame(load(File{format"CSV"}(string(exePath, "SoilTriangle.csv"))))
        end
        for i in 1:length(triangle.DP)
            if isBetween(sand, triangle.SandLow[i], triangle.SandHigh[i]) && isBetween(silt, triangle.SiltLow[i], triangle.SiltHigh[i]) && isBetween(clay, triangle.ClayLow[i], triangle.ClayHigh[i])
                return [triangle.Ksat[i], triangle.Sav[i], triangle.DP[i]]
            end
            #return error(string("Unable to classify soil texture"))
        end
    else
        return error(string("Soil textures are in percent, and must add to 100. Current sum is ", sum([sand, silt, clay])))
    end
end

# This takes the runoff from the PRZM .zts file and writes the .isd file for VFSMOD
function writeSediment(RUNF0, ESLS0, scen, fns)
    #ESLS0 comes in in tonnes/ha/day WRONG #############
    # ESL0 comes in tonnes/day
    solid = ESLS0 * constants().gInATonne #now in grams/day
    #scen.fieldAreaInHa*
    #RUNF0 comes in cm/ha/day
    #liquid = RUNF0/constants().cmInAMeter*scen.fieldAreaInHa*constants().mSqInAHa*constants().cm3InAM3#now in cm^3/day
    liquid = RUNF0 * scen.fieldAreaInHa * constants().mSqInAHa * constants().cmSqInAMSq #now in cm^3/day

    CI = solid / liquid #g/cm^3
    sediment = sedimentParameters(CI=CI)
    isd = open(fns.isdOutName, "w")
    write(isd, string("   ", sediment.NPART, "  ", @sprintf("%.4f", scen.COARSE), "  ", @sprintf("%.10f", sediment.CI), "  ", @sprintf("%.4f", scen.POR)), "\n")
    #Note that this currently writes sediment.SG (just the density of silicon) rather than scen.ρSediment
    write(isd, string("   ", @sprintf("%.4f", scen.DP), "  ", @sprintf("%.4f", scen.ρSed)))
    close(isd)
end

# This function takes a OPEN file (pontentially partially read), reads the next line,
# parses the first element of that line as a Float and returns it
function readfirst(F)
    x = split(readline(F), ',', limit=2)[1]
    x == "" ? 0.0 : parse(Float64, x)
end

# This function takes a OPEN file (pontentially partially read), reads the next line,
# parses the first element of that line as a Float and returns it
function readfirstSpace(F)
    x = split(readline(F), ' ', limit=10)
    while x[1] == ""
        popfirst!(x)
    end
    x == "" ? 0.0 : parse(Float64, x[1])
end

#The writes the .ikw file, and also returns filterStrip struct which gets passed to several other functions
function writeFilterStrip(width::Float64, scen::scenarioParameters, shapeFlag::Int64, prj::String, fn)
    # the way that VFSMOD conceptualizes the strip is 90° from PRZM
    # i.e., width and length are reversed

    # This calculation is based on one from the example sptreadsheet, with no citation or reason presented
    # nodes = width/3*11/2
    # As per the 2022 guidance
    nodes = width / 0.15

    # There is a minimum value of 11
    if nodes < 11
        nodes = 11
    end

    #Any values are rounded to the nearest odd number
    #Even integers are always increased
    nodes = floor(nodes / 2) * 2 + 1

    if shapeFlag == 1 # a square field with a rectangular pond and VFS
        # The length of the strip is just the square root of the field area - the hydraulic lentgh of the field is ignored
        length = sqrt(scen.fieldAreaInHa * constants().mSqInAHa)
    elseif shapeFlag == 2 # a round VFS around a round pond
        # Length calculated as the area over the width
        # Area is the difference between the circle within the outer edge ofthe strip and the area of the pond
        length = (pi * (sqrt(scen.pondAreaInM2 / pi) + width)^2 - scen.pondAreaInM2) / width
        # Old length was middle of strip
        #length = (sqrt(scen.pondAreaInM2/pi)+width/2)*2*pi
    elseif shapeFlag == 3 # A square VFS around a square pond
        # Length calculated as the area over the width
        # Area is the difference between the square within the outder edge of the strip and the area of the pond
        length = ((sqrt(scen.pondAreaInM2) + 2 * width)^2 - scen.pondAreaInM2) / width
    elseif shapeFlag == 4 # a squarish field dictates the length of the interface between the pond and a rectangular VFS on one side
        length = scen.fieldAreaInHa * constants().mSqInAHa / scen.hydraulicLengthInM
    else # this gives a square pond with a field on only one side, but this becomes the default
        # 
        length = sqrt(scen.pondAreaInM2)
    end

    filter = filterParameters(VL=width, FWIDTH=length, N=nodes, areaInMSq=width * length, areaInHa=width * length / constants().mSqInAHa, SX=width)
    ikw = open(fn, "w")
    write(ikw, string(prj, " - ", trunc(Int64, filter.VL), " m filter strip"), "\n") #, @sprintf("%.4f",scen.COARSE), "  ", @sprintf("%.10f",sediment.CI),"  ", @sprintf("%.4f",scen.POR)),"\n")
    write(ikw, string(@sprintf("%.3f", filter.FWIDTH)), "\n")
    write(ikw, string(@sprintf("%4s", filter.VL), @sprintf("%5s", filter.N), @sprintf("%5s", filter.θw), @sprintf("%5s", filter.CR), @sprintf("%6s", filter.MAXITER), @sprintf("%4s", filter.NPOL), @sprintf("%4s", filter.IELOUT), @sprintf("%4s", filter.KPG)), "\n")
    write(ikw, string(@sprintf("%2s", filter.NPROP)), "\n")
    write(ikw, string(@sprintf("%4s", filter.SX), @sprintf("%5s", filter.RNA), @sprintf("%5s", filter.SOA)), "\n")
    write(ikw, string(@sprintf("%2s", filter.IWQ)), "\n")
    close(ikw)
    return filter
end

#This writes the .igr file for VFSMOD
function writeGrass(fn)
    grass = grassParameters()
    igr = open(fn, "w")
    write(igr, string(" ", @sprintf("%.2f", grass.SS), @sprintf("%7s", grass.VN), @sprintf("%5s", grass.H), @sprintf("%6s", grass.VN2), @sprintf("%4s", grass.ICO)), "\n")
    write(igr, "--------------------------------------------", "\n")
    write(igr, " SS(cm)  Vn(s/cm^1/3)  H(cm)  Vn2(s/m^1/3)  ICO(0 or 1)", "\n")
    close(igr)
end

# This function takes the path to an SWI file and an existing
# waterQualityParameters struct and overwrites the
# chemical information in the struct with the values from the SWI
# It is called only once, before the main loop is run
function readChemicalParameters(fn)#, chemStruct::chemicalParameters)
    swi = open(fn, "r")
    #Throw away the first three lines
    for i in 1:3
        readline(swi)
    end
    KOCorKD = readline(swi) #4
    IKD = KOCorKD == "True" ? 1 : 2 # NEED TO VERIFY KD is 2, not 0 **********************************************************************
    KOCD = readfirst(swi) # 5 This line is a csv string of all KOCs or KDs, for parent, daughter, granddaughter
    for i = 6:12
        readline(swi) #Throw away seven more lines
    end
    refRate = readfirst(swi) #13
    refTemp = readfirst(swi) #14
    for i = 15:27
        readline(swi) #Throw away 13 more lines
    end
    Q10 = parse(Float64, readline(swi)) #28
    DGGHALF = refRate * Q10^((refTemp - constants().expectedRefTempInC) / 10)
    close(swi)
    chemStruct = chemicalParameters(IKD=IKD, KOCD=KOCD, DGGHALF=DGGHALF)
    return chemStruct
end

# This takes all sorts of stuff in order to write the .iwq file for VFSMOD
function writeWaterQualityParameters(inBetweenDays, RFLX1, EFLX1, thetaIn, chem::chemicalParameters, scen::scenarioParameters, water::waterQualityParameters, airTempIn, day, fns::inOutFileNames)
    #Only need one datum from this file, but it's position changes depending on the number of days of degradation after the event - that is also in the file

    owq = open(fns.owqFileName, "r")
    owqLine = ""
    owqWord = ""

    while owqWord != "area:" #scoot to right area (to the normalized data),
        #because there are two lines that end in "(mresn, IMOD=1)"
        owqLine = readline(owq)
        if length(split(owqLine)) > 0
            owqWord = split(owqLine)[end]
        end
    end

    while owqWord != "IMOB=1)" && owqWord != "days)" #scoot to right line, different by VFSMOD version
        owqLine = readline(owq)
        if length(split(owqLine)) > 0
            owqWord = split(owqLine)[end]
        end
    end

    mixingLayerResidue = parse(Float64, split(owqLine)[1]) # in mg/m2

    close(owq)

    massInLiquidInMg = RFLX1 * scen.fieldAreaInHa * constants().mSqInAHa * constants().cmSqInAMSq * constants().mgInAGram #starts as g/cm2/day
    massInSolidInMg = EFLX1 * scen.fieldAreaInHa * constants().mSqInAHa * constants().cmSqInAMSq * constants().mgInAGram #starts in g/cm2/day
    totalMassInMg = sum([massInLiquidInMg, massInSolidInMg])
    dgPin = totalMassInMg / (scen.fieldAreaInHa * constants().mSqInAHa)

    #In one scenario, it doesn't rain for 503 days, but VFSMOD will only let stuff degrade for 366 days
    if inBetweenDays > 365
        inBetweenDays = 365
    end

    airTemps = airTempIn[day:day+inBetweenDays-1]
    θs = thetaIn[day:day+inBetweenDays-1]

    iwq = open(fns.iwqOutName, "w")
    write(iwq, string(" ", water.IQPRO), "\n")
    write(iwq, string(" ", chem.IKD, " ", chem.KOCD, " ", scen.OCP,), "\n")
    write(iwq, string(" ", scen.CCP), "\n")
    write(iwq, string(" ", water.IDG), "\n")
    write(iwq, string(" ", inBetweenDays, " ", chem.DGGHALF, " ", scen.FC, " ", dgPin, " ", water.dgML, " ", chem.DGML, " ", mixingLayerResidue), "\n")
    for i in 1:length(airTemps)
        write(iwq, string(" ", airTemps[i]))
    end
    write(iwq, "\n")
    for i in 1:length(θs)
        write(iwq, string(" ", θs[i]))
    end
    write(iwq, "\n")
    write(iwq, string(" ", water.IMOB), "\n")
    close(iwq)
end

# The main loop looks for an OWQ file to get the precedent conditions. For the first iteration
# there isn't one, since VFSMOD hasn't been run yet.
function writeDummyOWQFile(fn)
    owq = open(fn, "w")
    for i in 1:12
        write(owq, "\n")
    end

    write(owq, string(0.0, " (mop)"), "\n")
    write(owq, string(0.0, " (mod)"), "\n")
    write(owq, string(0.0, " area:"), "\n")
    write(owq, string(0.0, " IMOB=1)"), "\n")

    close(owq)
end

# This takes the outputs from VFSMOD and  returns a waterQualityOutput struct that can be written directly to the .zts file
function readWaterQuality(scen::scenarioParameters, filter::filterParameters, fns::inOutFileNames)

    osm = open(fns.osmFileName, "r")
    owq = open(fns.owqFileName, "r")
    osmLine = ""

    while osmLine != "OUTPUTS" # scoot to the good stuff
        osmLine = readline(osm)
        if length(split(osmLine)) > 0
            osmLine = split(osmLine)[1]
        end
    end

    for i in 1:5
        readline(osm) #throw away 5 lines
    end

    runoffInM3 = parse(Float64, split(readline(osm))[5])
    for i in 1:11
        readline(osm) #throw away 11 lines
    end

    sedimentInG = parse(Float64, split(readline(osm))[6])
    owqLine = ""
    owqWord = ""

    # Ignore the entire file until the first line that ends in (mop)
    while owqWord != "(mop)"
        owqLine = readline(owq)
        if length(split(owqLine)) > 0
            owqWord = split(owqLine)[end]
        end
    end

    # That line has the pesticide mass on the solids leaving the strip
    pesticideSolidInmg = parse(Float64, split(owqLine)[1])
    # The next line has the pesticide mass in the dissolved phase leaving the strip
    pesticideWaterInmg = parse(Float64, split(readline(owq))[1])

    water = runoffInM3 * constants().cmInAMeter / filter.areaInMSq # now in cm/day
    solids = sedimentInG / constants().gInATonne # now in tonnes/day
    pesticideInWater = pesticideWaterInmg / constants().mgInAGram / (filter.areaInHa * constants().mSqInAHa * constants().cmSqInAMSq) #now in g/cm2/day
    pesticideInSolids = pesticideSolidInmg / constants().mgInAGram / (filter.areaInHa * constants().mSqInAHa * constants().cmSqInAMSq) #now in g/cm2/day

    close(osm)
    close(owq)
    # Write to a struct
    forZts = waterQualityOutput(RUNF0=water, ESLS0=solids, RFLX1=pesticideInWater, EFLX1=pesticideInSolids)
    return forZts
end

# This is a tool that only needs to be used once for each weather station
# It is used to generate initial soil moisture before each storm using PRZM
# It is not called anywhere in this code, and must be called manually from the REPL

# This newer version works with either PWC 1.52 (PRZM 5.02) or PWC 2.0x (PRZM 5.08)
# It only replaces what is necessary and doesn't require an external reference file
# thetaPath is the path to the folder containing the PWC generated PRZM5.inp file run on a scenario with the weather of interprets
# exePath is the path to the folder containing PRZM5.exe - if left blank, the copy looks in the thetaPath directory
function writePRZMTheta(thetaPath, exePath="", curveNumber=74, useDefaultRunoff=true)
    # Make a copy of the PRZM input file, since the output file will have the same name
    cp(string(thetaPath, "PRZM5.inp"), string(thetaPath, "PreTHETA.txt"), force=true)
    if !isfile(string(thetaPath, "PRZM5.inp"))
        error(string("Failed to locate file at ", thetaPath, "PRZM5.inp"))
        # return
    end
    oldInp = open(string(thetaPath, "PreTHETA.txt"), "r") #Open the copy of the PRZM5.inp as the source of most of the output
    VFSPM = open(string(thetaPath, "PRZM5.inp"), "w") #Open a new PRZM5.inp that will generate the VFS Soil Moisture from PRZM

    # Need to search through the original input file, line by line, using key "words"
    word = ""
    line = readline(oldInp) # First line of the old file
    write(VFSPM, line, "\n") # Write that down

    # Just copy and paste the original .inp file until the Record 1
    # Record 1 will be used to determine which version of PRZM (5.02 or 5.08) is being used
    while word != "1:" # The second word in "***Record 1:..." when split on spaces
        line = readline(oldInp) # Next line of the old file
        write(VFSPM, line, "\n") # Write that down
        if length(split(line, " ", keepempty=false)) > 1 # Want the second word, but only if the line has a space in it
            word = split(line, " ", keepempty=false)[2] # get the second word
        else
            word = ""
        end
    end

    if split(line, " ", keepempty=false)[end] == "anetd"
        PWC2 = true # In PWC 2, inicrp is not input in record 1, so the last word in this line is anetd
    else
        PWC2 = false
    end

    # PRZM 5.02 and 5.08 input files record the crop and runoff descriptors differently
    # In PRZM 5.02, it starts at Record 5
    if !PWC2
        while word != "5:"
            line = readline(oldInp) # Next line of the old file
            write(VFSPM, line, "\n") # Write that down
            if length(split(line, " ", keepempty=false)) > 1 # Want the second word, but only if the line has a space in it
                word = split(line, " ", keepempty=false)[2] # get the second word
            else
                word = ""
            end
        end

        line = "1, 0.1, 10,100, 0,   5" # Crop ID, canopyHoldup, rootDepth, canopyCover, WFMAX (optional), canopyHeight for turf
        write(VFSPM, line, "\n") # Write that down
        readline(oldInp) # Throw away the old crop information

        if useDefaultRunoff # Use the turf descriptors from the Turf_PA scenario

            # Write out the detailed runoff parameters from the Turf_PA scenario, all in one go
            line = string("***Record 6", "\n", "1, 26, 0", "\n", "***Record 7", "\n", "104,1604,105,1605,106,1506,1606,107,1607,108,1608,109,1609,110,1610,111,1611,112,1612,101,1601,102,1602,103,1503,1603,", "\n", "***Record 8", "\n", ".006,.002,.007,.004,.002,.007,.005,.003,.001,.005,.003,.003,.005,.009,.013,.013,.014,.014,.015,.015,.015,.015,.015,.015,.017,.012,", "\n", "***Record 9", "\n", ".110,.110,.110,.110,.110,.110,.110,.110,.110,.110,.110,.110,.110,.110,.110,.110,.110,.110,.110,.110,.110,.110,.110,.110,.110,.110,", "\n", "***Record 10", "\n", "74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,", "\n")
            write(VFSPM, line, "\n") # Write that down
            for i in 1:10
                readline(oldInp) # Throw away the old runoff information
            end
        else # Only overwrite the curve number
            # In PRZM 5.02, curve numbers are a list at Record 10
            while line != "***Record 10"
                line = readline(oldInp) # Next line of the old file
                write(VFSPM, line, "\n") # Write that down
            end

            curveNumber = string(curveNumber, ", ")
            curveNumbers = ""
            for i in 1:length(split(readline(oldInp), ",", keepempty=false)) # This line terminates in a comma, so -1 isn't necessary
                curveNumbers = string(curveNumbers, curveNumber)
            end
            write(VFSPM, string(curveNumbers), "\n")
        end

        line = readline(oldInp) # Get the header for the next record

    # In PRZM 5.08, curve number is the last value of a line or lines in Record 5
    else
        while word != "changes"
            line = readline(oldInp)
            write(VFSPM, line, "\n") # Write that down
            word = split(line, " ", keepempty=false)[end]
        end

        if useDefaultRunoff # use the runoff parameters from the Turf_PA

            line = "26, 0" # Number of hydro-event changes
            write(VFSPM, line, "\n") # Write that down
            readline(oldInp)
            line = readline(oldInp) # Header for the next record
            write(VFSPM, line, "\n") # Write that down

            # The set of hydr-event changes as a string
            line = string(" 1, 4,1972,0.006, 0.11, 74.0", "\n", "16, 4,    ,0.002, 0.11, 74.0", "\n", " 1, 5,    ,0.007, 0.11, 74.0", "\n", "16, 5,    ,0.004, 0.11, 74.0", "\n", " 1, 6,    ,0.002, 0.11, 74.0", "\n", "15, 6,    ,0.007, 0.11, 74.0", "\n", "16, 6,    ,0.005, 0.11, 74.0", "\n", " 1, 7,    ,0.003, 0.11, 74.0", "\n", "16, 7,    ,0.001, 0.11, 74.0", "\n", " 1, 8,    ,0.005, 0.11, 74.0", "\n", "16, 8,    ,0.003, 0.11, 74.0", "\n", " 1, 9,    ,0.003, 0.11, 74.0", "\n", "16, 9,    ,0.005, 0.11, 74.0", "\n", " 1,10,    ,0.009, 0.11, 74.0", "\n", "16,10,    ,0.013, 0.11, 74.0", "\n", " 1,11,    ,0.013, 0.11, 74.0", "\n", "16,11,    ,0.014, 0.11, 74.0", "\n", " 1,12,    ,0.014, 0.11, 74.0", "\n", "16,12,    ,0.015, 0.11, 74.0", "\n", " 1, 1,    ,0.015, 0.11, 74.0", "\n", "16, 1,    ,0.015, 0.11, 74.0", "\n", " 1, 2,    ,0.015, 0.11, 74.0", "\n", "16, 2,    ,0.015, 0.11, 74.0", "\n", " 1, 3,    ,0.015, 0.11, 74.0", "\n", "15, 3,    ,0.017, 0.11, 74.0", "\n", "16, 3,    ,0.012, 0.11, 74.0    ")
            write(VFSPM, line, "\n") # Write that down
            while word != "6:"
                line = readline(oldInp)
                word = split(line, " ", keepempty=false)[2]
            end
        else # Keep the everything but the curve number

            line = readline(oldInp) # Record 4
            write(VFSPM, line, "\n") # Write that down
            line = readline(oldInp) # Header for Record 5
            write(VFSPM, line, "\n") # Write that down

            line = readline(oldInp)
            word = split(line, " ", keepempty=false)[1]
            while word != "***Record"
                words = split(line, ",")
                line = string(words[1], ",", words[2], ",", words[3], ",", words[4], ",", words[5], ",", curveNumber)
                write(VFSPM, line, "\n") # Write that down
                line = readline(oldInp)
                word = split(line, " ", keepempty=false)[1]
            end
        end
    end

    # The number of crop periods is the next record
    write(VFSPM, line, "\n") # Write down the record header

    cropPeriods = parse(Int64, readline(oldInp)) # read the number of crop periods
    write(VFSPM, string(cropPeriods), "\n") # Write to the output file
    write(VFSPM, readline(oldInp), "\n") # Write whatever the record number and description is of the crop parameters

    # grab the first year of the crop descriptors
    firstYear = parse(Int64, split(readline(oldInp), ",")[3]) # Which is the 3rd element
    # ignore the rest of the crop parameters in the input file
    for i = 2:cropPeriods
        readline(oldInp)
    end

    # What to write for the crop parameters depends on the PWC version
    if PWC2
        lineEnd = ",      10,      0.6,        5,      0.1,        1"
    else
        lineEnd = ",       1"
    end

    # Write the same crop parameters for every year
    for i = 1:cropPeriods
        year = firstYear - 1 + i
        line = string(" 1, 4,", year, ",  15, 4,", year, ",   1,11,", year, lineEnd)
        write(VFSPM, line, "\n")
    end

    i = 0
    # Just copy and paste from the old .inp file until the last line of the file, but don't write that line
    while word != "Chemicals" # The first word of the last line of the section instructing PRZM what to put in the .zts file
        line = readline(oldInp) # Next line of the old file
        write(VFSPM, line, "\n") # Write that down
        word = split(line, " ")[end] # get the last word
    end

    line = PWC2 ? " 1, 1, False, 0, 0, 0, 0" : "       1       1" # Just one application and one chemical will do, but different for PRZM versions
    write(VFSPM, line, "\n")
    readline(oldInp) # skip the number of applications and chemicals in the old file

    write(VFSPM, readline(oldInp), "\n") #Copy the ***Record C2 line
    write(VFSPM, readline(oldInp), "\n") #Copy the first application

    # Throw away the rest of the applications in the old file
    while word != "***Record"
        line = readline(oldInp)
        word = split(line, " ", keepempty=false)[1]
    end

    # Just copy and paste from the old .inp file until the last line of the file, but don't write that line
    #while word != "INFL" # The first word of the last line of the section instructing PRZM what to put in the .zts file
    #    write(VFSPM, line, "\n") # Write that down
    #    line = readline(oldInp) # Next line of the old file
    #    word = split(line, ",")[1] # get the first word
    #end

    # Copy and paste until record U1, the number of lines of instructions for the output
    while word != "U1"
        line = readline(oldInp)
        write(VFSPM, line, "\n") # Write that down
        word = split(line, " ", keepempty=false)[end]
    end
    outputLines = parse(Int64, readline(oldInp))
    write(VFSPM, string(outputLines + 1), "\n") # Adding an extra line of output

    # copy to the end of the old file
    while word != "INFL" # The first word of the last line of the section instructing PRZM what to put in the .zts file
        line = readline(oldInp) # Next line of the old file
        write(VFSPM, line, "\n") # Write that down
        word = split(line, ",")[1] # get the first word
    end

    # Tell PRZM to write THETA as the last output column
    write(VFSPM, string("THET,0,TSER,   100,  100,    1.0"), "\n")

    # Close the files
    close(oldInp)
    close(VFSPM)

    # Run PRZM
    if exePath == ""
        turfPa = thetaPath[1:(end-1)]
        przm = `$turfPa\\PRZM5.exe`
        cd(thetaPath)
    else
        exePa = exePath[1:(end-1)]
        przm = `$exePa\\PRZM5.exe`
        cd(thetaPath)
    end
    run(przm)
end

# The old version combines elements of the 50 year crop run and a pre-made turf-specific run
# It works only with PWC 1.52's PRZM Input Files
function writePRZMTurf(turfPath)
    # Make a safe copy of the PRZM input file
    cp(string(turfPath, "PRZM5.inp"), string(turfPath, "PRZM5Crop.inp"), force=true)
    turfIn = open(string(turfPath, "PRZM5Turf.inp"), "r") # a file that holds the 
    crop = open(string(turfPath, "PRZM5Crop.inp"), "r")
    turf = open(string(turfPath, "PRZM5.inp"), "w")
    for i = 1:3
        write(turf, readline(crop), "\n")
        readline(turfIn)
    end
    line = readline(crop) #4
    readline(turfIn)
    #write(turf, string("Z:", line[26:end]),"\n") #4
    write(turf, line, "\n") #4
    write(turf, readline(crop), "\n") #5
    readline(turfIn)
    line = readline(crop) #6
    readline(turfIn)
    #write(turf, string("Z:", line[26:end]),"\n") #6
    write(turf, line, "\n") #6
    for i = 7:16
        write(turf, readline(crop), "\n")
        readline(turfIn)
    end
    for i = 17:81
        write(turf, readline(turfIn), "\n")
        readline(crop)
    end
    for i = 82:90
        write(turf, readline(crop), "\n")
    end
    horizons = trunc(Int64, readfirst(crop)) #91

    write(turf, string(horizons), "\n") #91
    for i = 92:(99+horizons)
        write(turf, readline(crop), "\n")
    end
    applicatations = trunc(Int64, readfirstSpace(crop)) #100 + horizons

    write(turf, string("     1     1"), "\n") #100 + horizons
    write(turf, readline(crop), "\n") #101 + horizons
    write(turf, readline(crop), "\n") #102 + horizons
    for i = 1:(applicatations-1)
        readline(crop) #just throwing these out
    end
    for i = 109:167
        write(turf, readline(crop), "\n")
    end
    # The last line of the file tells PRZM to write THET0 (soil moisture) instead of INFL0 into the .zts file
    write(turf, string("THET,0,TSER,   100,  100,    1.0"), "\n")
    close(crop)
    close(turf)
    close(turfIn)
    turfPa = turfPath[1:(end-1)]
    przm = `$turfPa\\PRZM5.exe`
    cd(turfPath)
    run(przm)
end

# This writes a new VVWMTransfer file by copying lines from the old one and inserting those that have changed
function writeVVWMTransfer(fileNames::inOutFileNames, filter::filterParameters, usInp::userInputs)
    vvwmIn = open(string(fileNames.vvwmTransferFileName[1:end-4], "_sansVFS.txt"), "r")
    vvwmOut = open(fileNames.vvwmTransferFileName, "w")
    readline(vvwmIn) #1
    write(vvwmOut, fileNames.przmOutName[1:end-4], "\n") #1
    for i = 2:29
        write(vvwmOut, readline(vvwmIn), "\n")
    end
    readline(vvwmIn)
    write(vvwmOut, string(fileNames.dailyWeatherFileName), "\n") #30
    for i = 31:58
        write(vvwmOut, readline(vvwmIn), "\n")
    end
    readline(vvwmIn) #59
    write(vvwmOut, string(trunc(Int64, filter.areaInMSq)), "\n") #59
    for i = 60:68
        write(vvwmOut, readline(vvwmIn), "\n")
    end
    seeEssVee = split(readline(vvwmIn), "\\")[end] #69, dude!
    oldCSV = string(workingPath, seeEssVee)
    newCSV = string(workingPath, Int64(usInp.stripWidthInM), "m_VFS_", seeEssVee)
    write(vvwmOut, string("\"", newCSV), "\n") #69, dude!
    for i = 70:71
        write(vvwmOut, string("\"", workingPath, "VFS_", split(readline(vvwmIn), "\\")[end]), "\n")
    end
    teeExTee = split(readline(vvwmIn), "\\")[end] #72
    oldTXT = string(workingPath, teeExTee)
    newTXT = string(workingPath, Int64(usInp.stripWidthInM), "m_VFS_", teeExTee)
    write(vvwmOut, string("\"", newTXT), "\n") #72
    for i = 73:83 #length(vvwm)
        write(vvwmOut, string("\"", workingPath, "VFS_", split(readline(vvwmIn), "\\")[end]), "\n")
    end
    if fileNames.isPWCFormat
        for i in 84:86
            write(vvwmOut, readline(vvwmIn), "\n")
        end
    end
    close(vvwmIn)
    close(vvwmOut)
    return [oldCSV[1:end-1], newCSV[1:end-1], oldTXT[1:end-1], newTXT[1:end-1]]
end

# This opens the output of VVWM and returns a dataframe that is easy to plot
# CURRENTLY NOT CALLED FROM VFSPipe
function getVVWMText(cn::columnNames, fn)
    vvwmTxt = open(fn, "r")

    # A little work-around for getting a DataFrame from a space-delimited format,using ReadCSV
    # Sadly, ReadCSV doesn't have a good skiplines, so a new, temporary space-delimited file without the unwanted header information is created first
    tempName = string(dirname(fn), "\\vtmp.csv")
    tempFile = open(tempName, "w")

    vvwmLine = "XXXX"
    #throw away the stuff that isn't wanted
    while vvwmLine[1:4] != "YEAR"
        vvwmLine = readline(vvwmTxt)
        if length(vvwmLine) < 4
            vvwmLine = "XXXX"
        end
    end
    #Get the first year of data before entering the loop
    vvwmLine = readline(vvwmTxt)
    #the line after the yearly EECs is a bunch of asterisks
    while vvwmLine[1:2] != "**"
        write(tempFile, vvwmLine, "\n")
        vvwmLine = readline(vvwmTxt)
    end
    close(vvwmTxt)
    close(tempFile)
    df = DataFrame(load(File{format"CSV"}(tempName), spacedelim=true, header_exists=false, colnames=cn.vvwmTxtColumns))
    rm(tempName)
    return df
end

# The bit that does everything
function vfsMain(usInp::userInputs)
    # Julia, or Juno, needs to have the working directory set
    # It also requires a drive change separate from the rest of the path
    cd(usInp.workingPath[1:2])
    cd(string(usInp.workingPath[3:end]))

    # All of the rest of the input and output file names can be generated programmatically
    inOutNames = writeFileNames(usInp)

    # Check if the expected files exist for an elegant exit with no dangling files
    if length(checkFilesExist(inOutNames, usInp)) > 0
        error(string("Unable to find: ", checkFilesExist(inOutNames, usInp)))
        # return
    end

    # The properties of the grass never change, but its existence is used for program control
    writeGrass(inOutNames.igrOutName)

    # Faster to write a dummy owq file than to test if it's the first run through the main loop every time
    writeDummyOWQFile(inOutNames.owqFileName)

    #Need a place to put the new .zts info
    przmOut = open(inOutNames.przmOutName, "w")
    #It contains the same header information. These lines are not conducive to a DataFrame, so just copy them as-is
    for i in 1:3
        write(przmOut, readlines(inOutNames.przmInName)[i], "\n")
    end

    #First, need the input data
    #Dataframes are easier with named columns
    cn = columnNames()
    #CSV.jl doesn't handle space delimited files well, so using CSVFiles

    # Check for the existence of the input files before trying to load them
    if !isfile(inOutNames.przmInName)
        error("Unable to locate ", przmInName)
        #return
    end
    przmIn = DataFrame(load(File{format"CSV"}(inOutNames.przmInName), spacedelim=true, skiplines_begin=3, header_exists=false, colnames=cn.ztsColumns))

    if !isfile(inOutNames.thetaInName)
        error(string("Unable to locate soil moisture file at ", inOutNames.thetaInName))
        #return
    else
        thetaIn = open(inOutNames.thetaInName)
        cNames = split(readlines(thetaIn)[3], " ", keepempty=false)
        close(thetaIn)
        if "THET0" in cNames
            thetaPosition = findall(x -> x == "THET0", cNames)[1]
        else
            error("Unable to find THET0 column in the soil moisture file ", inOutFileNames.thetaInName)
            #return
        end
    end

    thetaIn = DataFrame(load(File{format"CSV"}(inOutNames.thetaInName), spacedelim=true, skiplines_begin=3, header_exists=false))[!, thetaPosition]

    #The Hourly Precipitation File (HPF), if being used
    if usInp.useHPF
        precipIn = DataFrame(load(File{format"CSV"}(inOutNames.hpfInName), spacedelim=true, header_exists=false, colnames=cn.hpfColumns))
    else
        # Only need the daily precipitation, but want it in a dataframe with a header
        # leading to this slightly awkward work-around
        if inOutNames.dailyWeatherFileName[end-2:end] == "wea"
            weatherIn = DataFrame(load(File{format"CSV"}(inOutNames.dailyWeatherFileName), spacedelim=false, header_exists=false, colnames=cn.weaColumns))
        else
            # Need to know how many columns the dvf file has, since it's not the same between PMRA and EPA
            dvfile = open(inOutNames.dailyWeatherFileName)
            dvfColumnCount = length(split(readline(dvfile), " ", keepempty = false))
            close(dvfile)
            if dvfColumnCount < 10
                weatherIn = DataFrame(load(File{format"CSV"}(inOutNames.dailyWeatherFileName), spacedelim=true, header_exists=false, colnames=cn.dvfColumns))
            else
                weatherIn = DataFrame(load(File{format"CSV"}(inOutNames.dailyWeatherFileName), spacedelim=true, header_exists=false, colnames=cn.dvf14Columns))
            end
        end
        precipIn = weatherIn[!, [:Total]]
        stormLength = usInp.stormLengthInHours * constants().secondsInAnHour
    end

    airTempIn = weatherIn[!, :Temperature]
    # Free up the memory
    weatherIn = []
    #Need the scenario information to be avaialable to several functions
    scenario = readScenarioParameters(inOutNames, usInp)
    chem = readChemicalParameters(inOutNames.swiFileName)

    # Properties of the filter strip do not change during the run
    filterStrip = writeFilterStrip(usInp.stripWidthInM, scenario, usInp.shapeFlag, "VFSM", inOutNames.ikwOutName)

    # Some parameters of the water do not change, and are assigned their own struct
    water = waterQualityParameters(IQPRO=usInp.pesticideEquation)

    #An index for keeping track of days bewtween simulations, to account for degradation
    inBetweenDays = 1

    # In order to keep the load in the VFS current, VFSMOD is run for every runoff event
    for day in 1:size(przmIn, 1)

        # Although VVWM probably doesn't read the dates - PRZM calls them 'dummy fields'
        # It's nice to keep the data and the format "just in case"
        yr = string(Int64(przmIn.Yr[day]))
        mo = Int8(przmIn.Mo[day])
        if mo < 10
            mo = string(" ", string(mo))
        else
            mo = string(mo)

        end
        dy = Int8(przmIn.Dy[day])
        if dy < 10
            dy = string(" ", string(dy))
        else
            dy = string(dy)
        end

        # If there is runoff, do all the stuff to run VFSMOD;
        #if there isn't, just copy the line from the old .zts file
        #Here we go!
        if przmIn.RUNF0[day] > 0
            #Need to know the number of days until the NEXT runoff event to be simulated
            #This is used by VFSMOD to calculate degradation in the strip between runoff events
            # The 2022 guidance clarifies that it INCLUDES the day of a runoff event, so it starts at 1
            inBetweenDays = 1
            # Need to prevent an end-of-read error if there is runoff on the last day
            if length(przmIn.RUNF0) - day > 0
                while przmIn.RUNF0[day+inBetweenDays] == 0
                    # Need to prevent an end-of-read error after the last runoff event too
                    inBetweenDays += 1
                    if inBetweenDays > length(przmIn.RUNF0) - day
                        break
                    end
                end
            end

            #Start with the rain - write a new irn file
            #The runoff calculation needs the total event time, either returned from the rain function, if using hourly weather,
            #Or set by the user if using daily weather information
            if usInp.useHPF
                eventTime, irn = writePrecipitation(precipIn[day, :][5:end], inOutNames)
            else
                eventTime = stormLength
                irn = writePrecipitation(precipIn.Total[day], stormLength, inOutNames)
            end
            #Next, the runoff from the field
            # 2022 Guidance indicates that runoff from rain gives an 8 hour hydrograph, but runoff from snowmelt should give a 12 hour hydrograph
            if precipIn.Total[day] > 0
                eventTime = stormLength
            else
                eventTime = 12
            end
            iro = writeRunoff(przmIn.RUNF0[day], scenario, filterStrip, eventTime, inOutNames)

            #write the soil file
            if day > 1
                writeSoil(thetaIn[(day-1)], scenario, inOutNames) #Want the water content from the day before
            else
                writeSoil(thetaIn[(day)], scenario, inOutNames) #The water content from the first day will have to do
            end

            #And the sediment file
            writeSediment(przmIn.RUNF0[day], przmIn.ESLS0[day], scenario, inOutNames)

            #And finally, the water quality
            writeWaterQualityParameters(inBetweenDays, przmIn.RFLX1[day], przmIn.EFLX1[day], thetaIn, chem, scenario, water, airTempIn, day, inOutNames)

            #Run VFSMOD
            run(inOutNames.vfsmod)
            println(string("Simulated Day: ", day))
            vfsOut = readWaterQuality(scenario, filterStrip, inOutNames)
            write(przmOut, string(yr, " ", mo, " ", dy, "         ", likePrzm(vfsOut.RUNF0), "   ", likePrzm(vfsOut.ESLS0), "   ", likePrzm(vfsOut.RFLX1), "   ", likePrzm(vfsOut.EFLX1), "   ", likePrzm(vfsOut.DCON1), "   ", likePrzm(vfsOut.INFL0)), "\n")

        else #There was no runoff, so no change in the .zts file
            write(przmOut, string(yr, " ", mo, " ", dy, "         ", likePrzm(przmIn.RUNF0[day]), "   ", likePrzm(przmIn.ESLS0[day]), "   ", likePrzm(przmIn.RFLX1[day]), "   ", likePrzm(przmIn.EFLX1[day]), "   ", likePrzm(przmIn.DCON1[day]), "   ", likePrzm(przmIn.INFL0[day])), "\n")

            #A day with no runoff is a day for degradation
            inBetweenDays += 1
        end
    end
    close(przmOut)

    if usInp.isLastRun
        rm(inOutNames.igrOutName)
    end
    # Writes the new VVWM Transfer File, and returns file names for plotting
    oldCSV, newCSV, oldTXT, newTXT = writeVVWMTransfer(inOutNames, filterStrip, usInp)

    return inOutNames.vvwm, oldTXT, newTXT
end
