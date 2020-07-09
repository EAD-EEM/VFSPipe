### This file contains all of the structs and the functions that manipulate them
# that are needed to run VFSMod from PRZM outputs, and VVWM from VFSMOD outputs

using CSVFiles, DataFrames, PrettyTables, Printf, CSV, Parameters, Plots #Gtk removed

#constants
@with_kw struct constants @deftype Int64
    #Not sure what's up
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

@with_kw struct soilParameters{R <:Real}
    #Using the variable names from VFSMOD
    θI::R = 0.0 #This is the only value that changes from run to run
    SM::R = 0.0 # Maximum surface storage
    SCHK::R = 0.5 # relative distance from upper edge..?
end

@with_kw struct chemicalParameters{R <:Real}
    IKD::Int8 = 1 #1=KOC used
    KOCD::R = 100.0 #
    DGGHALF::R = 30.0 #Half=life, not K
end
#These are the parameters required by the water quality file, and are not sourced from other structs
# Other structs required for the output include chemical and scenario parameters
@with_kw struct waterQualityParameters{R <:Real}
    IQPRO::Int8 = 1 #1 = Based on Sabbagh
    IDG::Int8 = 2 #2=degradation by decay rate only
    dgML::R = 2.0 #Surface mixing layer thickness = constant
end

@with_kw struct sedimentParameters{R <:Real}
    #Using the variable names from VFSMOD
    NPART::Int8 = 7 # 7 = read user input values
    CI::R = 0.0 #This is the only value value changes from run to run
    SG::R = 2.65 #This is the density of quartz
end

@with_kw struct scenarioParameters{R <:Real}
    #For a few calculations
    fieldAreaInHa::R = 10.0 #in ha
    pondAreaInM2::R = 10000.0 #in m^2
    #For water quality
    OCP::R = 3.2 #Organic Carbon %
    CCP::R = 19.0 #Clay %
    FC::R = 0.31 # Field capacity
    # For sediments
    COARSE::R = 0.75 # really just sand fraction
    POR::R = 0.472 # same as saturated water content
    DP::R = 0.0098 #This should be modified for each SCENARIO
    #For Soil
    VKS::R = 6.06E-6 # Saturated hydraulic conductivity (should be calculated from scenario)
    SAV::R = 0.1101 # Green-Ampts' suction at wet front (should be calculated from scenario)
    #Note in appedix says 'equivent to POR value'
    θSoil::R = 0.47 # Saturated soil-water content (should be taken from scenario, same as porosity for now)
    ρSed::R = 2.65 # Density of sediment particles
    # Slope needed for filterParameters
    SOA::R = 0.0200 # (should be read)
end
@with_kw struct waterQualityOutput @deftype Float64
    RUNF0 = 0.0 #cm per ha per day - water per field area
    ESLS0 = 0.0 #tonnes per ha per day - solids per field area
    RFLX1 = 0.0 # g per ha per day - mass of pesticide by runoff
    EFLX1 = 0.0 # g per ha per day - mass of pesticide by erosion
    DCON1 = 0.0 # g per ha per day - mass of degradate 1 by runoff
    INFL0 = 0.0 # g per ha per day - mass of degradate 1 by erosion
end

@with_kw struct filterParameters @deftype Float64
    FWIDTH = 363.821 # Width of strinp in m
    VL = 3.0 # Length of strip in m
    N::Int64 = 11  # Number of nodes
    θw = 0.5 # Time weighted factor, Crank-Nicholson solution
    CR = 0.8 # Courant number
    MAXITER::Int64 = 350 # Maximum number of interations
    NPOL::Int64 = 3 # Number of nodal points
    IELOUT::Int8 = 1 #Flag for output
    KPG::Int8 = 1 # Flag for Petrox-Galerkin solution
    NPROP::Int8 = 1 #Number of segments
    SX = 3.0 #Distance of filter of uniform surface properties - default is same as VL
    RNA = 0.45 #Manning's n for segment
    SOA = 0.0200 # Slope not read here - write method takes from scenarioParameters
    IWQ::Int8 = 1 #Flag for water quality calculation, 1 = do it
    areaInMSq = 1091.5 # Area in m^2
    areaInHa = 0.1091 # Area in Ha
end

@with_kw struct grassParameters @deftype Float64
    SS = 2.15
    VN = 0.012
    H = 18
    VN2 = 0.05
    ICO::Int8 = 0
end

@with_kw struct inOutFileNames @deftype String

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

@with_kw struct columnNames @deftype Array{String}
    ztsColumns = ["Yr","Mo","Dy","RUNF0","ESLS0","RFLX1","EFLX1","DCON1","INFL0"]
    hpfColumns = ["Date","Duration","Total","Midnight","OneAM","TwoAM","ThreeAM","FourAM","FiveAM","SixAM","SevenAM","EightAM","NineAM","TenAM","ElevenAM","Noon","OnePM","TwoPM","ThreePM","FourPM","FivePM","SixPM","SevenPM","EightPM","NinePM","TenPM","ElevenPM"]
    dvfColumns = ["Date", "Total", "A", "Temperature", "B", "C"]
    vvwmTxtColumns = ["Run", "Peak", "FourDay","TwentyOneDay","SixtyDay","NinetyDay","OneYear","PWPeak","PWTwentyOneDay"]
    vvwmCsvColumns = ["Depth","Average","PWAverage","Peak"]
end

#The .zts file uses an inconvenient format for the numbers, required for the correct spacing
#This function just formats numbers that way - in base 10 scientific notation with one digit left of the
#decimal, 4 digits right of the decimal, and three digits in the exponent
function likePrzm(num)
    return string(@sprintf("%.4E", num)[1:8], @sprintf("%03d", parse(Int64, @sprintf("%.4E", num)[9:10])))
end

#Write all the file names
function writeFileNames(workingPath, exePath, turfPath, projectName, pwcName, twoCharacterCode)
    if length(projectName) > 6
        error("Project name must be <= 6 characters")
    end
    executeString = string(projectName,twoCharacterCode,".prj")
    exePath = exePath[1:end-1] #To make the path 'command friedly'
    vfsmod = `$exePath\\vfsm.exe $executeString`
    vvwmTransferFileName = string(workingPath, "vvwmTransfer.txt")
    vvwm = `$exePath\\vvwm.exe $vvwmTransferFileName`
    swiFileName = string(workingPath, pwcName, ".SWI")
    if isfile(swiFileName)
        # Do nothing
    else
        swiFileName = string(workingPath, pwcName, ".PWC")
    end
    przmInName = string(workingPath, projectName, "MN.zts")
    #przmInName = string(workingPath, projectName, "NM.zts")
    cp(string(workingPath,pwcName,".zts"), przmInName, force = true)
    przmOutName = string(workingPath,projectName, ".zts")
    #For now, since we don't really need an HPF - could create from dvf or met right here
    #In future, should explore real hourly with the 'new' met data
    hpfInName = string(workingPath, projectName, ".HPF")
    prjInName = string(workingPath, projectName, ".prj")
    prjOutName = string(workingPath, projectName, twoCharacterCode, ".prj")
    scnFileName = string(workingPath, readlines(swiFileName)[49], ".SCN")
    dailyWeatherFileName = string(workingPath, split(readlines(swiFileName)[50], "\\")[end])
    thetaInName = string(turfPath,split(readlines(swiFileName)[50], "\\")[end][1:(end-4)], "Turf.zts")

    #Create a project file for VFSMOD to read the file names too
    fileTypes = ["ikw","iso","igr","isd","irn","iro","iwq","og1","og2","ohy","osm","osp","owq"]
    prj = open(prjOutName, "w")
    for i in 1:length(fileTypes)
        write(prj, string(fileTypes[i],"=",projectName,twoCharacterCode,".",fileTypes[i]),"\n")
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
    owqFileName = string(workingPath, readlines(prjOutName)[13][5:end])
    osmFileName = string(workingPath, readlines(prjOutName)[11][5:end])

    #Place all of these in a struct so they cannot be accidentally overwritten
    inOutNames = inOutFileNames(vfsmod = vfsmod, vvwm = vvwm, przmInName = przmInName, thetaInName = thetaInName, przmOutName = przmOutName, hpfInName = hpfInName, prjInName = prjInName, prjOutName = prjOutName, swiFileName = swiFileName, scnFileName = scnFileName, dailyWeatherFileName = dailyWeatherFileName, ikwOutName = ikwOutName, irnOutName = irnOutName, iroOutName = iroOutName, isoOutName = isoOutName, isdOutName = isdOutName, iwqOutName = iwqOutName, owqFileName = owqFileName, igrOutName = igrOutName, osmFileName = osmFileName, vvwmTransferFileName = vvwmTransferFileName)
    return inOutNames
end

#This takes a vector of HOURLY rain intensities and writes out an irn file of durations and intensities
#This version trims zeroes off both ends to reduce run time of VFSMOD
function writePrecipitation(row::DataFrame,fn)
    #Containers
    times = Vector{Int16}()
    rates = Vector{Float32}()
    #Convert this row from a DataFrame into a vector for easy manipulation
    row = vec(permutedims(Vector(row)))
    #Remove any zeroes from the beginning
    while row[1] == 0
        deleteat!(row,1)
    end
    #Remove any zeroes from the end
    while row[end] == 0
        deleteat!(row,length(row))
    end

    #Do not want reduduncancy, but easier to write them all and remove than to control what to write
    for i in 1:length(row)
        t = (i-1)*constants().secondsInAnHour # We want the hour, not the index
        times = [times;t]
        r = row[i]/constants().cmInAMeter/constants().secondsInAnHour #HPF is in cm/h, irn is in m/s
        rates = [rates;r]
    end
    i = 1 #length(times)
    while length(times[i:end]) > 2
        if rates[i] == rates[i+1] #&& rates[i] == rates[i+2]
            deleteat!(times, i+1)
            deleteat!(rates, i+1)
        else
            i += 1
        end
    end

    #The runoff calculation needs the precipitation time
    #Need to record it before adding on the zeroes
    stormTime = maximum(times)

    #According to the UFSMOD manual, it is generally it is not raining at the end of the simulation
    #The example pipe appears to end the simulation with a two-hour rain-less period
    if maximum(times) < (constants().hoursInADay-2)*constants().secondsInAnHour
        times = [times; (maximum(times)+1)] #one second after the last non-zero
        rates = [rates;0]
        times = [times;((maximum(times)-1)+2*constants().secondsInAnHour)] #two hours after the last non-zero
        rates = [rates;0]
    end

    #VFSMOD expects a description of the data at the top of the data
    pushfirst!(times, length(times))
    pushfirst!(rates, maximum(rates))
    irn = open(fn, "w")
    for j in 1:length(times)
        write(irn, string(@sprintf("%9s",times[j]), "  ", @sprintf("%.4E",rates[j])), "\n")
        #write(irn, string(" ", @sprintf("%08u",times[j]), "  ", @sprintf("%.4E",rates[j])), "\n")
    end
    close(irn)
    return stormTime
end

# This version takes a daily rain amount and generates an stormTime seconds long storm from it
function writePrecipitation(precip::Float64, stormTime::Int64, fn)

    numberOfDataInAFakeStorm = 4
    startOfStorm = 0
    oneSecondAfterStorm = stormTime +1
    secondlyRate = precip/constants().cmInAMeter/stormTime #.dvf is in cm/storm, rainfall intensity in .irn is in m/s
    # THere are two hours of simulation after the storm
    endOfSimulationTime = stormTime + 2*constants().secondsInAnHour
    # Unless, of course, the storm lasts for > 22 hours
    if endOfSimulationTime > constants().hoursInADay*constants().secondsInAnHour
        endOfSimulationTime = constants().hoursInADay*constants().secondsInAnHour
    end

    irn = open(fn, "w")
    # The first line has the number of events to read and the maximum rate, which is the same as the rate, since the hydrograph is square
    write(irn, string(@sprintf("%9s",numberOfDataInAFakeStorm), "  ", @sprintf("%.4E",secondlyRate)), "\n")
    write(irn, string(@sprintf("%9s",startOfStorm), "  ", @sprintf("%.4E",secondlyRate)), "\n")
    write(irn, string(@sprintf("%9s",stormTime), "  ", @sprintf("%.4E",secondlyRate)), "\n")
    write(irn, string(@sprintf("%9s",oneSecondAfterStorm), "  ", @sprintf("%.4E",0.0)), "\n")
    write(irn, string(@sprintf("%9s",endOfSimulationTime), "  ", @sprintf("%.4E",0.0)), "\n")
    close(irn)
end

#This version just transmutes 23 hours of precipitation rates into an irn file
#rather, it will if I finish it
# function getRainTime(row)
#     times = Vector{Int16}()
#     rates = Vector{Float32}()
#     for hour in 1:length(row) #We want the hour, not the index
#         t = (hour-1)*constants().secondsInAnHour
#         times = [times;t]
#         r = row[hour]
#         rates = [rates;r]
#     end
#     df = DataFrame()
#     df.Times = times
#     df.Rates = rates
#
#     #According to the UFSMOD manual, it is generally it is not raining at the end of the simulation
#     #Simulation ends at midnight to write the ZTS line
#     if max(df.Times) < constants()hoursInADay*constants().secondsInAnHour - 1
#         push!(df, ((constants()hoursInADay*constants().secondsInAnHour-1), 0))
#     end
#     return df
# end

function writeSoil(θ, scen::scenarioParameters, fn)
    soil = soilParameters(θI = θ)
    iso = open(fn, "w")
    write(iso, string("  ", @sprintf("%.2E",scen.VKS),"  ", @sprintf("%.4f",scen.SAV),"  " ,@sprintf("%.4f",scen.θSoil),"  " ,@sprintf("%.4f",soil.θI),"  ", @sprintf("%.4f",soil.SM),"  ", @sprintf("%.4f",soil.SCHK)))
    close(iso)
end
#This function takes the daily runoff from PRZM, a storm length, and a filemane
#and writes an iro file
#There are many assumptions built in, brought forward from the example pipe
function writeRunoff(RUNF0,scen::scenarioParameters,filter::filterParameters, eventLength,fn)

    fieldAreaInMSq = scen.fieldAreaInHa*constants().mSqInAHa #in ha - confirmed in VVWM manual that runoff is normalized to FIELD area (not waterbody area)

    #Current defaults - here, required to back-calculate the runoff from the field which VFSMOD will then use to calculate in m/s
    #SWIDTH = 363.821 #in meters`
    #SLENGTH = 274.861 #in meters

    SWIDTH = filter.FWIDTH#(sqrt(filter.areaInMSq/π)+filter.FWIDTH)*2π
    SLENGTH = fieldAreaInMSq/SWIDTH

    timeStep::Int64 = (constants().secondsInAnHour/2) # Arbitrary choice - half an hour
    #The 8 hour storms in the example have 10 time steps - it is not known why
    #For fun, here, the time step is set and the number derivied
    numberOfSteps::Int16 = eventLength/timeStep
    #Get from cm/ha/day to m^3/day, and then per second
    runOn = RUNF0/constants().cmInAMeter*fieldAreaInMSq/eventLength
    #Containers
    times = Vector{Int16}()
    rates = Vector{Float32}()

    #write the times and the rate for the event
    for i in 0:numberOfSteps
        push!(times,i*timeStep)
        push!(rates,runOn)
    end

    #The second line has the number of timesteps and the peak rate of the hydrograph
    pushfirst!(times, numberOfSteps)
    pushfirst!(rates, maximum(rates))

    iro = open(fn, "w")
    write(iro, string(" ",SWIDTH," ",SLENGTH), "\n")
    for i = 1:length(times)
        write(iro, string(@sprintf("%9s",times[i]), " ", @sprintf("%.4E", rates[i])), "\n")
    end
    close(iro)
end

function readScenarioParameters(fn)
    scn = open(fn, "r")
    #Throw away the first three lines
    for i in 1:6
        readline(scn)
    end
    fieldArea = 10.0
    fldArea = readline(scn) #7
    if fldArea != ""
        fieldArea = parse(Float64, fldArea)/constants().mSqInAHa
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
    slope = readfirst(scn)/constants().cent #49 - must be fractional, not percent
    for i in 50:52
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
    POR = 1-(ρ/constants().ρQuartz)
    texture = getSoilClass(sandPercent,(100-sum([sandPercent,clayPercent])), clayPercent)
    #This derivation probably highlights a logical flaw in determining the bulk density of sediments
    #In PWC scenario development rather than providing a better value for sediment bulk density
    #It shall be left here but not used in the writing of the sediment input file for now
    ρSedimentParticle = (ρBenthicBulk - benthicPorosity*constants().ρWater)/(1 - benthicPorosity)
    scenStruct = scenarioParameters(fieldAreaInHa = fieldArea, OCP = ocPercent, CCP = clayPercent, FC = fieldCapacity, COARSE = sandPercent/constants().cent, POR = POR, VKS = texture[1], SAV = texture[2], DP = texture[3], θSoil = POR, SOA = slope) #ρSed = ρSedimentParticle,
    return scenStruct
end

# Determines if b is between a and c
#Type not declared, so it works on strings too
function isBetween(b,a,c)
    if a <= b <=c
        return true
    else
        return false
    end
end

#This function takes percentages of sand, silt, and clay in that order, determines the class
#And returns a vector of Ksat, Sav, and DP in that order
function getSoilClass(sand::Real,silt::Real,clay::Real)
    # Data from http://www.fao.org/tempref/FI/CDrom/FAO_Training/FAO_Training/General/x6706e/x6706e06.htm
    # Actually, that data has an error - for Silty Loam the range of Silt was wrong
    # table 4
    #Check that the values add to 100, but with some leeway for rounding error
    if isBetween(sum([sand,silt,clay]), 99.5, 100.5)
        triangle = DataFrame(load(File(format"CSV", string(exePath, "SoilTriangle.csv"))))
        for i in 1:length(triangle.DP)
            if isBetween(sand, triangle.SandLow[i], triangle.SandHigh[i]) && isBetween(silt, triangle.SiltLow[i], triangle.SiltHigh[i]) && isBetween(clay, triangle.ClayLow[i], triangle.ClayHigh[i])
                return [triangle.Ksat[i], triangle.Sav[i],triangle.DP[i]]
            end
        #return error(string("Unable to classify soil texture"))
        end
    else
        return error(string("Soil textures are in percent, and must add to 100. Current sum is ", sum([sand,silt,clay])))
    end
end

function writeSediment(RUNF0,ESLS0, scen, fn)
    #ESLS0 comes in in tonnes/ha/day WRONG #############
    # ESL0 comes in tonnes/day
    solid = ESLS0*constants().gInATonne #now in grams/day
    #scen.fieldAreaInHa*
    #RUNF0 comes in cm/ha/day
    #liquid = RUNF0/constants().cmInAMeter*scen.fieldAreaInHa*constants().mSqInAHa*constants().cm3InAM3#now in cm^3/day
    liquid = RUNF0*scen.fieldAreaInHa*constants().mSqInAHa*constants().cmSqInAMSq #now in cm^3/day

    CI = solid/liquid #g/cm^3
    sediment = sedimentParameters(CI = CI)
    isd = open(fn, "w")
    write(isd, string("   ", sediment.NPART, "  ", @sprintf("%.4f",scen.COARSE), "  ", @sprintf("%.10f",sediment.CI),"  ", @sprintf("%.4f",scen.POR)),"\n")
    #Note that this currently writes sediment.SG (just the density of silicon) rather than scen.ρSediment
    write(isd, string("   ", @sprintf("%.4f",scen.DP),"  ", @sprintf("%.4f",scen.ρSed)))
    close(isd)
end

# This function takes a OPEN file (pontentially partially read), reads the next line,
# parses the first element of that line as a Float and returns it
function readfirst(F)
  x = split(readline(F),',', limit=2)[1]
  x == "" ? 0.0 : parse(Float64, x)
end

# This function takes a OPEN file (pontentially partially read), reads the next line,
# parses the first element of that line as a Float and returns it
function readfirstSpace(F)
  x = split(readline(F),' ', limit=10)
  while x[1] == ""
      popfirst!(x)
  end
  x == "" ? 0.0 : parse(Float64, x[1])
end

function writeFilterStrip(width::Float64, scen::scenarioParameters, prj::String, fn)
    # the way that VFSMOD conceptualizes the strip is 90° from PRZM
    # i.e., width and length are reversed

    # This calculation is based on one from the example sptreadsheet, with no citation or reason presented
    # It may appear in the VFSMOD guidance
    nodes = width/3*11/2

    # There appears to be a minimum value of 11
    # This method is slightly biased in that it rounds all even integers up
    if nodes < 11
        nodes = 11
    end

    #Any values are rounded to the nearest odd number
    nodes = floor(nodes/2)*2+1

    length = (sqrt(scen.pondAreaInM2/pi)+width/2)*2*pi
    filter = filterParameters(VL = width, FWIDTH = length, N = nodes, areaInMSq = width*length, areaInHa = width*length/constants().mSqInAHa, SX = width)
    ikw = open(fn, "w")
    write(ikw, string(prj," - ", trunc(Int64, filter.VL), " m filter strip"),"\n") #, @sprintf("%.4f",scen.COARSE), "  ", @sprintf("%.10f",sediment.CI),"  ", @sprintf("%.4f",scen.POR)),"\n")
    write(ikw, string(@sprintf("%.3f", filter.FWIDTH)),"\n")
    write(ikw, string(@sprintf("%4s", filter.VL), @sprintf("%5s", filter.N), @sprintf("%5s", filter.θw), @sprintf("%5s", filter.CR), @sprintf("%6s", filter.MAXITER), @sprintf("%4s", filter.NPOL), @sprintf("%4s", filter.IELOUT), @sprintf("%4s", filter.KPG)), "\n")
    write(ikw, string(@sprintf("%2s", filter.NPROP)),"\n")
    write(ikw, string(@sprintf("%4s", filter.SX), @sprintf("%5s", filter.RNA), @sprintf("%5s", filter.SOA)),"\n")
    write(ikw, string(@sprintf("%2s", filter.IWQ)),"\n")
    close(ikw)
    return filter
end

function writeGrass(fn)
    grass = grassParameters()
    igr = open(fn, "w")
    write(igr, string(" ", @sprintf("%.2f", grass.SS), @sprintf("%7s", grass.VN), @sprintf("%5s", grass.H), @sprintf("%6s", grass.VN2), @sprintf("%4s", grass.ICO)),"\n")
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
    DGGHALF = refRate*Q10^((refTemp-constants().expectedRefTempInC)/10)
    close(swi)
    chemStruct = chemicalParameters(IKD = IKD, KOCD = KOCD, DGGHALF = DGGHALF)
    return chemStruct
end

function writeWaterQualityParameters(inBetweenDays, RFLX1, EFLX1, thetaIn, chem::chemicalParameters, scen::scenarioParameters, airTempIn, day, owqfn, fn)
    #Only need one datum from this file, but it's position changes depending on the number of days of degradation after the event - that is also in the file
    # oldInBetweenDays = parse(Int64,split(readlines(owqfn)[13])[6])

    # If there are more than 100, VFSMOD adds an extra blank line
    # if oldInBetweenDays > 100
    #     oldInBetweenDays += 1
    # end

    # There are 37 other lines not required that do not change
    # massInFilterRead = split(readlines(owqfn)[37+oldInBetweenDays])[1]
    # if massInFilterRead[1] == "-"[1]
    #     massInFilter = 0.0
    # else
    #     massInFilter = parse(Float64,massInFilterRead) # in mg
    # end
    # if massInFilter < 0
    #     println(massInFilter, RFLX1, EFLX1)
    #     massInFilter = 0.0 # in mg

    # end
    # massInFilterRead = split(readlines(owqfn)[37+oldInBetweenDays])[1]
    # if length(massInFilterRead) > 4
    #     if massInFilterRead[end-3] == "+"[1]
    #         massInFilter = parse(Float64, string(massInFilterRead[1:end-4],"E",massInFilterRead[end-3:end]))
    #     else
    #         massInFilter = parse(Float64,massInFilterRead) # in mg
    #     end
    # else
    #     massInFilter = parse(Float64,massInFilterRead) # in mg
    # end
    owq = open(owqfn, "r")
    owqLine = ""

    while owqLine != "Outputs"
        owqLine = readline(owq)
        if length(split(owqLine)) > 0
            owqLine = split(owqLine)[1]
        end
    end

    for i in 1:19
        readline(owq) #throw away 5 lines
    end

    # pesticideWaterInmg = parse(Float64, split(readlines(owqfn)[52+oldInBetweenDays])[1])
    # pesticideSolidInmg = parse(Float64, split(readlines(owqfn)[51+oldInBetweenDays])[1])

    massInFilter = parse(Float64, split(readline(owq))[1])
    close(owq)

    if isnan(massInFilter)
        massInFilter = 0.0
    end
    #massInLiquidInMg = RFLX1*scen.fieldAreaInHa*constants().mgInAGram #starts as g/ha/day
    massInLiquidInMg = RFLX1*scen.fieldAreaInHa*constants().mSqInAHa*constants().cmSqInAMSq*constants().mgInAGram #starts as g/cm2/day
    #massInSolidInMg = EFLX1*scen.fieldAreaInHa*constants().mgInAGram #starts in g/ha/day
    #massInSolidInMg = EFLX1*constants().mgInAGram #starts in g/day
    massInSolidInMg = EFLX1*scen.fieldAreaInHa*constants().mSqInAHa*constants().cmSqInAMSq*constants().mgInAGram #starts in g/cm2/day
    totalMassInMg = sum([massInFilter,massInLiquidInMg,massInSolidInMg])
    dgPin = totalMassInMg/(scen.fieldAreaInHa*constants().mSqInAHa)
    #dgPin = totalMassInMg/(363.9155481418726*3)
    #For Debugging only
    #inBetweenDays = 0

    #In one scenario, it doesn't rain for 503 days, but VFSMOD will only let stuff degrade for 366 days
    if inBetweenDays > 365
        inBetweenDays = 365
    end

    airTemps = airTempIn[day:day+inBetweenDays]
    θs = thetaIn[day:day+inBetweenDays]
    #water = waterQualityParameters(dgPin = totalMass)
    water = waterQualityParameters()
    iwq = open(fn,"w")
    write(iwq, string(" ", water.IQPRO),"\n")
    write(iwq, string(" ", chem.IKD, " ", chem.KOCD, " ", scen.OCP, ), "\n")
    write(iwq, string(" ", scen.CCP), "\n")
    write(iwq, string(" ", water.IDG), "\n")
    write(iwq, string(" ", inBetweenDays, " ", chem.DGGHALF, " ", scen.FC, " ", dgPin, " ", water.dgML), "\n")
    for i in 1:length(airTemps)
        write(iwq, string(" ", airTemps[i]))
    end
    write(iwq, "\n")
    for i in 1:length(θs)
        write(iwq, string(" ", θs[i]))
    end
    write(iwq, "\n")
    close(iwq)
end

# The main loop looks for an OWQ file to get the precedent conditions. For the first iteration
# there isn't one, since VFSMOD hasn't been run yet.
function writeDummyOWQFile(fn)
    owq = open(fn, "w")
    for i in 1:12
        write(owq, "\n")
    end
    write(owq, "    No. of days between events=    0", "\n") #13

    for i in 14:16
        write(owq, "\n")
    end
    write(owq, " Outputs for Water Quality", "\n") #17

    for i in 18:36
        write(owq, "\n")
    end
    write(owq, "            0.0 mg = Pesticide surface residue at next event (after degradation,  0 days)") #37


    close(owq)
end

function readWaterQuality(scen::scenarioParameters, filter::filterParameters, owqfn, osmfn)
    #oldInBetweenDays = parse(Int64,split(readlines(owqfn)[13])[6])

    osm = open(osmfn, "r")
    owq = open(owqfn, "r")
    # If there are more than 100, VFSMOD adds an extra blank line
    #if oldInBetweenDays > 100
        #oldInBetweenDays += 1
    #end

    # Sometimes VFSMOD doesn't write the .og1 file so the OSM file is missing 20 lines
    # This is probably a problem
    # diffCheck = readlines(osmfn)[86]
    # if diffCheck == " OUTPUTS"
    #     missingOgFile = 20
    # else
    #     missingOgFile = 0
    # end

    osmLine = ""

    while osmLine != "OUTPUTS"
        osmLine = readline(osm)
        if length(split(osmLine)) > 0
            osmLine = split(osmLine)[1]
        end
    end

    for i in 1:5
        readline(osm) #throw away 5 lines
    end

    runoffInM3 = parse(Float64, split(readline(osm))[5])
    #runoffInM3 = parse(Float64, split(readlines(osmfn)[112-missingOgFile])[5])
    #sedimentinKg = parse(Float64, split(readlines(owqfn)[20+oldInBetweenDays])[1])*(1 - parse(Float64, split(readlines(owqfn)[23+oldInBetweenDays])[1])/constants().cent) # percent to fraction
    #println(string("Runoff is: ", runoffInM3,"m^3"))
    for i in 1:11
        readline(osm) #throw away 11 lines
    end

    sedimentInG = parse(Float64, split(readline(osm))[6])
    #println(string("Sediment mass is: ", sedimentInG,"g"))
    #sedimentRead = split(readlines(osmfn)[124-missingOgFile])[6]
    # If there's crazy erosion, outflow is calculated as infinity
    #Hopefully this was only required during debugging
    # if sedimentRead == "Infinity"
    #     sedimentInG = parse(Float64,split(readlines(osmfn)[123-missingOgFile])[6])
    # elseif sedimentRead[end-3] == "+"[1]
    #     sedimentInG = parse(Float64, string(sedimentRead[1:end-4],"E",sedimentRead[end-3:end]))
    # else
    #     sedimentInG = parse(Float64,sedimentRead)#split(readlines(osmfn)[124])[6])
    # end
    #
    owqLine = ""

    while owqLine != "Outputs"
        owqLine = readline(owq)
        if length(split(owqLine)) > 0
            owqLine = split(owqLine)[1]
        end
    end

    for i in 1:33
        readline(owq) #throw away 5 lines
    end

    # pesticideWaterInmg = parse(Float64, split(readlines(owqfn)[52+oldInBetweenDays])[1])
    # pesticideSolidInmg = parse(Float64, split(readlines(owqfn)[51+oldInBetweenDays])[1])

    pesticideWaterInmg = parse(Float64, split(readline(owq))[1])
    pesticideSolidInmg = parse(Float64, split(readline(owq))[1])

    if isnan(runoffInM3)
        runoffInM3 = 0.0
    end
    if isnan(sedimentInG)
        sedimentInG = 0.0
    end
    if isnan(pesticideWaterInmg)
        pesticideWaterInmg = 0.0
    end
    if isnan(pesticideSolidInmg)
        pesticideSolidInmg = 0.0
    end

    # if isnan(parameter)
    #     for parameter in [runoffinM3,sedimentinKg,pesticideWaterInmg,pesticideSolidInmg]
    #         parameter = 0.0
    #     end
    # end

    #Convert to the correct dimensions and units
    # water = runoffInM3*constants().cmInAMeter/scen.fieldAreaInHa/constants().mSqInAHa # now in cm/ha/day
    # solids = sedimentInG/scen.fieldAreaInHa/constants().gInATonne # now in tonnes/ha/day
    # pesticideInWater = pesticideWaterInmg/constants().mgInAGram/scen.fieldAreaInHa #now in g/ha/day
    # pesticideInSolids = pesticideSolidInmg/constants().mgInAGram/scen.fieldAreaInHa #now in g/ha/day

    # To match the example, the area of the FILTER is used for the output to zts
    # That is matched when the VVWM transfer file is written.
    # VVWM just turns the per area quantities back into raw volumes and masses, so it doesn't matter, as long as they match
    # The filter strip area is in m^2
    water = runoffInM3*constants().cmInAMeter/filter.areaInMSq # now in cm/day
    #solids = sedimentInG/constants().gInATonne/filter.areaInHa # now in tonnes/ha/day
    solids = sedimentInG/constants().gInATonne # now in tonnes/day
    #pesticideInWater = pesticideWaterInmg/constants().mgInAGram/filter.areaInHa #now in g/ha/day
    pesticideInWater = pesticideWaterInmg/constants().mgInAGram/(filter.areaInHa*constants().mSqInAHa*constants().cmSqInAMSq) #now in g/cm2/day
    #pesticideInSolids = pesticideSolidInmg/constants().mgInAGram/filter.areaInHa #now in g/ha/day
    pesticideInSolids = pesticideSolidInmg/constants().mgInAGram/(filter.areaInHa*constants().mSqInAHa*constants().cmSqInAMSq) #now in g/cm2/day

    close(osm)
    close(owq)
    # Write to a struct
    forZts = waterQualityOutput(RUNF0 = water, ESLS0 = solids, RFLX1 = pesticideInWater, EFLX1 = pesticideInSolids)
    return forZts
end

# This is a tool that only needs to be used once on each turf scenario
# It is not colled anywhere in this code
function writePRZMTurf(turfPath) #Do you need these?
    # Make a safe copy - just in case?
    cp(string(turfPath, "PRZM5.inp"), string(turfPath,"PRZM5Crop.inp"), force = true)
    turfIn = open(string(turfPath, "PRZM5Turf.inp"), "r")
    crop = open(string(turfPath,"PRZM5Crop.inp"), "r")
    turf = open(string(turfPath, "PRZM5.inp"), "w")
    for i = 1:3
        write(turf, readline(crop), "\n")
        readline(turfIn)
    end
    line = readline(crop) #4
    readline(turfIn)
    write(turf, string("Z:", line[26:end]),"\n") #4
    write(turf, readline(crop), "\n") #5
    readline(turfIn)
    line = readline(crop) #6
    readline(turfIn)
    write(turf, string("Z:", line[26:end]),"\n") #6
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
    println(horizons)
    write(turf, string(horizons), "\n") #91
    for i = 92:(99+horizons)
        write(turf, readline(crop), "\n")
    end
    applicatations = trunc(Int64, readfirstSpace(crop)) #100 + horizons
    println(applicatations)
    write(turf,string("     1     1"), "\n" ) #100 + horizons
    write(turf, readline(crop), "\n") #101 + horizons
    write(turf, readline(crop), "\n") #102 + horizons
    for i = 1:(applicatations-1)
        readline(crop) #just throwing these out
    end
    for i = 109:169
        write(turf, readline(crop), "\n")
    end
    write(turf, string("THET,0,TSER,   100,  100,    1.0"))
    close(crop)
    close(turf)
    close(turfIn)
    turfPa = turfPath[1:(end-1)]
    przm = `$turfPa\\PRZM5.exe`
    cd(turfPath)
    run(przm)
end

# function writeVVWMTransfer(fileNames::inOutFileNames, filter::filterParameters, runYears)
#     # Make a safe copy - just in case?
#     cp(string(fileNames.vvwmTransferFileName), string(fileNames.vvwmTransferFileName[1:end-4], "sansVFS.txt"), force = true)
#     vvwmIn = open(string(fileNames.vvwmTransferFileName[1:end-4], "sansVFS.txt"), "r")
#     vvwmOut = open(fileNames.vvwmTransferFileName, "w")
#     readline(vvwmIn) #1
#     write(vvwmOut, fileNames.przmOutName[1:end-4], "\n") #1
#     for i = 2:29
#         write(vvwmOut, readline(vvwmIn), "\n")
#     end
#     readline(vvwmIn)
#     write(vvwmOut, string(fileNames.dailyWeatherFileName), "\n") #30
#     for i = 31:55
#         write(vvwmOut, readline(vvwmIn), "\n")
#     end
#     readline(vvwmIn) #56
#     write(vvwmOut, string(runYears), "\n") #56
#     incomingIn = split(readline(vvwmIn), ",")#57
#     incomingOut = Vector{Int64}()
#     for i in 1:runYears
#         push!(incomingOut, parse(Int64, incomingIn[i]))
#     end
#     write(vvwmOut, join(map(string,incomingOut),","), "\n")#57
#     #write(vvwmOut, readline(vvwmIn), "\n")#57
#     write(vvwmOut, readline(vvwmIn), "\n")#58
#     readline(vvwmIn) #59
#     write(vvwmOut, string(trunc(Int64, filter.areaInMSq)), "\n") #59
#     for i = 60:68
#         write(vvwmOut, readline(vvwmIn), "\n")
#     end
#     for i = 69:83 #length(vvwm)
#         write(vvwmOut, string("\"",workingPath,split(readline(vvwmIn), "\\")[end]), "\n")
#     end
#     close(vvwmIn)
#     close(vvwmOut)
# end

function writeVVWMTransfer(fileNames::inOutFileNames, filter::filterParameters)
    # Make a safe copy - allows direct comparison of what enters the pond
    cp(string(fileNames.vvwmTransferFileName), string(fileNames.vvwmTransferFileName[1:end-4], "sansVFS.txt"), force = true)
    vvwmIn = open(string(fileNames.vvwmTransferFileName[1:end-4], "sansVFS.txt"), "r")
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
    newCSV = string(workingPath, "VFS_", seeEssVee)
    write(vvwmOut, string("\"", newCSV), "\n") #69. dude!
    for i = 70:71
        write(vvwmOut, string("\"",workingPath,"VFS_",split(readline(vvwmIn), "\\")[end]), "\n")
    end
    teeExTee = split(readline(vvwmIn), "\\")[end] #72
    oldTXT = string(workingPath, teeExTee)
    newTXT = string(workingPath, "VFS_", teeExTee)
    write(vvwmOut, string("\"", newTXT), "\n") #72
    for i = 73:83 #length(vvwm)
        write(vvwmOut, string("\"",workingPath,"VFS_",split(readline(vvwmIn), "\\")[end]), "\n")
    end
    close(vvwmIn)
    close(vvwmOut)
    return [oldCSV[1:end-1],newCSV[1:end-1],oldTXT[1:end-1],newTXT[1:end-1]]
end

## This function is no longer called
## It was from debugging, when a local copy was required
# function runVVWM(fileNames::inOutFileNames)
#     cp(fileNames.vvwm, string(workingPath, "vvwmTemp.exe"), force = true)
#     workingPathTrim = string(workingPath)[1:end-1]
#     vvwmExe = `$workingPathTrim\\vvwmTemp.exe`
#     run(vvwmExe)
#     rm(string(workingPath, "vvwmTemp.exe"))
# end

function getVVWMText(cn::columnNames, fn)
    vvwmTxt = open(fn, "r")

    # A little work-around for getting a DataFrame from a space-delimited format,using ReadCSV
    # Sadly, ReadCSV doesn't have a good skiplines, so a new, temporary space-delimited file without the unwanted header information is created first
    tempName = string(dirname(fn),"\\vtmp.csv")
    tempFile = open(tempName, "w")

    # The OLD way was not elegant and didn't handle different numbers of years of run
    # for i in 1:19
    #     readline(vvwmTxt)
    # end
    # for i in 1:simYears
    #     write(tempFile, readline(vvwmTxt), "\n")
    # end

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
    df = DataFrame(load(File(format"CSV", tempName), spacedelim = true, header_exists = false, colnames = cn.vvwmTxtColumns))
    rm(tempName)
    return df
    #return DataFrame(load(File(format"CSV", tempName), spacedelim = true, header_exists = false, colnames = cn.vvwmTxtColumns))
end

@with_kw struct userInputs @deftype String
    projectName = ""
    stripWidthInM::Float64 = 0.0
    twoCharacterCode = ""
    workingPath = ""
    pwcName = ""
    useHPF::Bool = false
    stormLengthInHours::Int64 = 8
    exePath = ""
    turfPath = ""
end

#Can save some work by following a rigid filename convention
# projectName = "SOILS" #Six characters only
# stripWidthInM = 10.0
# #Can replace the last two letters of the input files for the temporary working files so they do not get written over
# twoCharacterCode = "RG"
# workingPath = "Z:\\SharedwithVM\\VFS\\Soils\\9005\\"
# pwcName = "Persistane-REG"
# # If there's a proper hourly precipitation file, use that, otherwise generate the precipitation event programmatically from the daily value
# useHPF = false
# stormLengthInHours = 8 #only read when useHPF is false
#
# #Less often changed are the paths to unchanging files
# exePath = "Z:\\SharedwithVM\\VFS\\executables\\"
# turfPath = "Z:\\SharedwithVM\\VFS\\CanadianTurfZts\\"

function vfsMain(usInp::userInputs)
    # Julia, or Juno, needs to have the working directory set
    # It also requires a drive change separate from the rest of the path
    #cd("Z:")
    cd(usInp.workingPath[1:2])
    cd(string(usInp.workingPath[3:end]))

    # All of the rest of the input and output file names can be generated programmatically
    inOutNames = writeFileNames(usInp.workingPath, usInp.exePath, usInp.turfPath, usInp.projectName, usInp.pwcName, usInp.twoCharacterCode)

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
    przmIn = DataFrame(load(File(format"CSV", inOutNames.przmInName), spacedelim = true, skiplines_begin = 3, header_exists = false, colnames = cn.ztsColumns))
    thetaIn = DataFrame(load(File(format"CSV", inOutNames.thetaInName), spacedelim = true, skiplines_begin = 3, header_exists = false))[!,9]
    #The Hourly Precipitation File (HPF), if being used
    if usInp.useHPF
        precipIn = DataFrame(load(File(format"CSV", inOutNames.hpfInName), spacedelim = true, header_exists = false, colnames = cn.hpfColumns))
    else
        # Only need the daily precipitation, but want it in a dataframe with a header
        # leading to this slightly awkward work-around
        precipIn = DataFrame(load(File(format"CSV", inOutNames.dailyWeatherFileName), spacedelim = true, header_exists = false, colnames = cn.dvfColumns))
        precipIn = precipIn[!, [:Total]]
        stormLength = usInp.stormLengthInHours*constants().secondsInAnHour
    end

    airTempIn = DataFrame(load(File(format"CSV", inOutNames.dailyWeatherFileName), spacedelim = true, header_exists = false))[!,4]
    #Need the scenario information to be avaialable to several methods
    scenario = readScenarioParameters(inOutNames.scnFileName)
    chem = readChemicalParameters(inOutNames.swiFileName)

    # Properties of the filter strip do not change during the run
    filterStrip = writeFilterStrip(usInp.stripWidthInM, scenario, usInp.projectName, inOutNames.ikwOutName)
    # The properties of the grass never change
    writeGrass(inOutNames.igrOutName)

    #An index for keeping track of days bewtween simulations, to account for degradation
    inBetweenDays = 0

    ### Turns out that THETO is the soil moisture at the BEGINNING of the simuation
    #Want the soil moisture from BEFORE the rain - but in the rare instance that it rains on the first day of the simulation, the first value must suffice. This variable is set to one at the bottom of the main for loop
    #turfOffset = 0

    # In order to keep the load in the VFS current, VFSMOD is run for every rain event in the HPF file
        #No longer True # Most of the time, there isn't runoff.
        #No longer True # Only when there is runoff does the VFS do anything
        # But precipitation PROBABLY causes infiltration of the pesticide...
    for day in 1:size(precipIn, 1)

        # Although VVWM probably doesn't read the dates - PRZM calls them 'dummy fields'
        # It's nice to keep the data and the format "just in case"
        yr = string(Int8(przmIn.Yr[day]))
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

        # If it rains, do all the stuff to run VFSMOD;
        #if it doesn't, just copy the line from the old .zts file
    #    if precipIn.Total[day] > 0 #Here we go!
        if przmIn.RUNF0[day] > 0
            #Need to know the number of days until the NEXT rain event to be simulated
            #This used by VFSMOD to calculate degradation in the strip between rain events
            inBetweenDays = 0
            # Need to prevent an end-of-read error if it rains on the last day
            if length(przmIn.RUNF0) - day > 0
                while przmIn.RUNF0[day+1+inBetweenDays] == 0
                    # Need to prevent an end-of=read error after the last rain event
                    inBetweenDays += 1
                    if inBetweenDays+1 > length(przmIn.RUNF0)-day
                        break
                    end
                end
            end

            #Start with the rain - write a new irn file
            #The runoff calculation needs the total event time, either returned from the rain function, if using hourly weather,
            #Or set by the user if using daily weather information
            if usInp.useHPF
                eventTime = writePrecipitation(precipIn[day,:][5:end], inOutNames.irnOutName)
            else
                eventTime = stormLength
                writePrecipitation(precipIn.Total[day], stormLength, inOutNames.irnOutName)
            end
            #Next, the runoff from the field
            writeRunoff(przmIn.RUNF0[day], scenario, filterStrip, eventTime, inOutNames.iroOutName)

            #write the soil file
            writeSoil(thetaIn[(day)], scenario, inOutNames.isoOutName)#Want the

            #And the sediment file
            writeSediment(przmIn.RUNF0[day], przmIn.ESLS0[day], scenario, inOutNames.isdOutName)

            #And finally, the water quality
            writeWaterQualityParameters(inBetweenDays, przmIn.RFLX1[day], przmIn.EFLX1[day], thetaIn, chem, scenario, airTempIn, day, inOutNames.owqFileName, inOutNames.iwqOutName)
            #Write a line to the new zts file
            #write(przmOut, "Hey...you should have run VFSMOD!", "\n")
            run(inOutNames.vfsmod)
            println(string("Simulation: ",day))
            vfsOut = readWaterQuality(scenario, filterStrip, inOutNames.owqFileName, inOutNames.osmFileName)
            write(przmOut, string(yr, " ", mo, " " , dy, "         ", likePrzm(vfsOut.RUNF0), "   ", likePrzm(vfsOut.ESLS0), "   ", likePrzm(vfsOut.RFLX1), "   ", likePrzm(vfsOut.EFLX1), "   ", likePrzm(vfsOut.DCON1), "   ", likePrzm(vfsOut.INFL0)), "\n")

        else #There was no rain, so no change in the VFS
            write(przmOut, string(yr, " ", mo, " " , dy, "         ", likePrzm(przmIn.RUNF0[day]), "   ", likePrzm(przmIn.ESLS0[day]), "   ", likePrzm(przmIn.RFLX1[day]), "   ", likePrzm(przmIn.EFLX1[day]), "   ", likePrzm(przmIn.DCON1[day]), "   ", likePrzm(przmIn.INFL0[day])), "\n")
            #write(przmOut, join(@sprintf("%.4E3",permutedims(Vector(przmIn[day,:])))," "), "\n")

            #A day with no precip is a day for degradation
            inBetweenDays += 1
        end
        ### Turns out that THETO is the soil moisture at the BEGINNING of the simuation
        # Again, this is only used to prevent a read-error, by trying to get the
        # soil moisture in the turf on the day before the przm simulation starts
        #global turfOffset = 1
    end
    close(przmOut)

    # Writes the new VVWM Transfer File, and returns file names for plotting
    oldCSV,newCSV,oldTXT,newTXT = writeVVWMTransfer(inOutNames, filterStrip)
    #writeVVWMTransfer(inOutNames, filterStrip,14)
    #run(inOutNames.vvwm)
    return inOutNames.vvwm, oldTXT, newTXT
end
# These last bits allow quick plotting of results
