#Generating the headers for a BDF file
#Author: Danick Lamoureux
#LM2 project under Frédérick Gosselin's supervision
#Date: 2022-05-05

#All headers are modified from Simcenter's original headers

def headers_hydro(file, P1, nmodes = 10):
    #Headers to use for full hydroelastic analysis. Requires modal and coupled acoustics analysis first


    file.write(r"""$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*                    Aeroelasticity approximation of hydrodynamic damping
$*                    for Simcenter Nastran version 2020.2
$*
$*
$ EXECUTIVE CONTROL
$*
ID, Aeroelastic analysis, added mass
SOL 145
COMPILE FLUTTER NOLIST NOREF $
$MALTER '$ COUNTER FOR NUMBER OF M/K PAIRS' $
ALTER 1 $""")
    for i in range(nmodes):
        file.write("\nMATGEN ,/PPP"+str(i+1)+"/1/1 $\n")
        file.write("ADD PPP"+str(i+1)+",/PP"+str(i+1)+"/"+str(P1[i,i])+" $")
        if i+1>1:
            file.write("\nMATGEN ,/ZERO"+str(i)+"/7/"+str(i)+"/1/0 $\n")
            file.write("TRNSP ZERO"+str(i)+"/ZEROT"+str(i)+" $\n")
            file.write("MATGEN ,/CP"+str(i)+"/6/"+str(i+1)+"/"+str(i)+"/1 $\n")
            if i+1==2:
                file.write("MERGE PP1,ZEROT1,ZERO1,PP2,CP1,/A1 $")
            else:
                file.write("MERGE A"+str(i-1)+",ZEROT"+str(i)+",ZERO"+str(i)+",PP"+str(i+1)+",CP"+str(i)+",/A"+str(i)+" $")
    file.write("""
EQUIVX A"""+str(nmodes-1)+"""/BETAMAT/ALWAYS $
$Mass matrix scaling
MPYAD MHH,BETAMAT,/MHHX $
MPYAD MHHX,BETAMAT,/MHHNEW $
EQUIVX MHHNEW/MHH/ALWAYS $
$Damping matrix scaling
MPYAD BHH,BETAMAT,/BHHX $
EQUIVX BHHX/BHH/ALWAYS $
ENDALTER
CEND
$*
$ CASE CONTROL
$*
ECHO = NONE
SPC = 100
METHOD = 101
FMETHOD = 200
CMETHOD = 102
OUTPUT
AEROF = ALL
APRESSURE = ALL
DISPLACEMENT(PLOT,REAL) = ALL
EDE(PLOT,AVERAGE,RMAG) = ALL
FORCE(PLOT,REAL,CENTER) = ALL
SPCFORCES(PLOT,REAL) = ALL
STRESS(PLOT,REAL,VONMISES,CENTER) = ALL
$*  Step: Subcase - Eigenvalue Method 1
SUBCASE 1
  LABEL = Subcase - Eigenvalue Method 1
  METHOD = 101
  CMETHOD = 102
$*
$* BULK DATA
$*
BEGIN BULK
$*
$* SOLUTION CARDS
$*
$*  Modeling Object: Real Eigenvalue - Lanczos1
EIGRL, 101, 10.0000, , """+str(nmodes)+""", 0, 7, , MASS
$* Modeling Object: Complex Eigenvalue - Lanczos
EIGC,102, CLAN, , , , , """+str(nmodes)+"""
$*
$* PARAM CARDS
$*
PARAM      K6ROT100.0000
PARAM     OIBULK     YES
PARAM    OMACHPR     YES
PARAM       POST      -2
PARAM    POSTEXT     YES
PARAM    UNITSYS   MN-MM
PARAM    BAILOUT      -1
PARAM, AUTOSPC, NO
$PARAM, EXTOUT, DMIGPCH
$*"""+'\n')

def headers_aero(file):
    #Headers to use for flutter analysis
    file.write(r"""$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*                    Aeroelasticity analysis of a hydrofoil
$*                    for Simcenter Nastran version 2020.2
$*
$*
$ EXECUTIVE CONTROL
$*
ID, Aeroelastic analysis, no added mass
SOL 145
CEND
$*
$* CASE CONTROL
$*
ECHO = NONE
SPC = 100
METHOD = 101
FMETHOD = 200
OUTPUT
AEROF = ALL
APRESSURE = ALL
DISPLACEMENT(PLOT,REAL) = ALL
EDE(PLOT,AVERAGE,RMAG) = ALL
FORCE(PLOT,REAL,CENTER) = ALL
SPCFORCES(PLOT,REAL) = ALL
STRESS(PLOT,REAL,VONMISES,CENTER) = ALL
$*  Step: Subcase - Eigenvalue Method 1
SUBCASE 1
  LABEL = Subcase - Eigenvalue Method 1
  METHOD = 101
$*
$* BULK DATA
$*
BEGIN BULK
$*
$* SOLUTION CARDS
$*
$*  Modeling Object: Real Eigenvalue - Lanczos1
EIGRL, 101, 10.0000, 10.0E3, 20, 0, 7, , MASS
$*
$* PARAM CARDS
$*
PARAM      K6ROT100.0000
PARAM     OIBULK     YES
PARAM    OMACHPR     YES
PARAM       POST      -2
PARAM    POSTEXT     YES
PARAM    UNITSYS   MN-MM
$*"""+'\n')

def headers_acoustics_real(file):
    #Headers to use for uncoupled acoustics analysis (obsolete)
    file.write(r"""$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*                    Acoustics analysis of a hydrofoil (uncoupled)
$*                    for Simcenter Nastran version 2020.2
$*
$*
$ EXECUTIVE CONTROL
$*
ID, Acoustics analysis, no added damping, uncoupled fluid-structure
SOL 103
CEND
$*
$* CASE CONTROL
$*
ECHO = NONE
SPC = 100
MODSEL = ALL
MODSEL(FLUID) = ALL
FLSTCNT ACSYM=NO ACOUT=PEAK ASCOUP=YES PREFDB=0.0 SKINOUT=FREEFACE ,
AGGPCH=YES SFEF70=NO
METHOD(STRUCTURE) = 101
METHOD(FLUID) = 101
OUTPUT
ACVELOCITY(PRINT,REAL) = ALL
BCRESULTS(TRACTION,FORCE,PLOT) = ALL
DISPLACEMENT(PLOT,REAL) = ALL
EDE(PLOT,AVERAGE,RMAG) = ALL
FORCE(PLOT,REAL,CENTER) = ALL
MEFFMASS(PRINT,SUMMARY) = YES
PRESSURE(PLOT,REAL,TOTAL) = ALL
SPCFORCES(PLOT,REAL) = ALL
STRESS(PLOT,REAL,VONMISES,CENTER) = ALL
$*  Step: Subcase - Eigenvalue Method 1
SUBCASE 1
  LABEL = Subcase - Eigenvalue Method 1
  METHOD(STRUCTURE) = 101
  METHOD(FLUID) = 101
  SPC = 100
  OUTPUT
  ACVELOCITY(PRINT,REAL) = ALL
  BCRESULTS(TRACTION,FORCE,PLOT) = ALL
  DISPLACEMENT(PLOT,REAL) = ALL
  EDE(PLOT,AVERAGE,RMAG) = ALL
  FORCE(PLOT,REAL,CENTER) = ALL
  MEFFMASS(PRINT,SUMMARY) = YES
  PRESSURE(PLOT,REAL,TOTAL) = ALL
  SPCFORCES(PLOT,REAL) = ALL
  STRESS(PLOT,REAL,VONMISES,CENTER) = ALL
$*
$* BULK DATA
$*
BEGIN BULK
$*
$* SOLUTION CARDS
$*
TABDMP1      300       G                                                +
+         0.00000.0000001.000+100.020000    ENDT
$*  Modeling Object: Real Eigenvalue - Lanczos
EIGRL, 101, 10.0000, 10.0E3, 20, 0, 7, , MASS
$*
$* PARAM CARDS
$*
PARAM      K6ROT100.0000
PARAM     OIBULK     YES
PARAM    OMACHPR     YES
PARAM       POST      -2
PARAM    POSTEXT     YES
PARAM    UNITSYS   MN-MM
$*"""+'\n')

def headers_acoustics_complex(file):
    #Headers to use for coupled acoustics analysis
    file.write(r"""$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*                    Acoustics analysis of a hydrofoil (coupled)
$*                    for Simcenter Nastran version 2020.2
$*
$*
$ EXECUTIVE CONTROL
$*
ID, Acoustics analysis, no added damping, coupled fluid-structure
SOL 107
CEND
$*
$* CASE CONTROL
$*
ECHO = NONE
SPC = 100
MODSEL = ALL
MODSEL(FLUID) = ALL
FLSTCNT ACSYM=NO ACOUT=PEAK ASCOUP=YES PREFDB=0.0 SKINOUT=FREEFACE ,
AGGPCH=YES SFEF70=NO
METHOD(STRUCTURE) = 101
METHOD(FLUID) = 101
CMETHOD = 102
OUTPUT
ACVELOCITY(PRINT,REAL) = ALL
BCRESULTS(TRACTION,FORCE,PLOT) = ALL
DISPLACEMENT(PLOT,REAL) = ALL
EDE(PLOT,AVERAGE,RMAG) = ALL
FORCE(PLOT,REAL,CENTER) = ALL
MEFFMASS(PRINT,SUMMARY) = YES
PRESSURE(PLOT,REAL,TOTAL) = ALL
SPCFORCES(PLOT,REAL) = ALL
STRESS(PLOT,REAL,VONMISES,CENTER) = ALL
$*  Step: Subcase - Eigenvalue Method 1
SUBCASE 1
  LABEL = Subcase - Eigenvalue Method 1
  METHOD(STRUCTURE) = 101
  METHOD(FLUID) = 101
  SPC = 100
  OUTPUT
  ACVELOCITY(PRINT,REAL) = ALL
  BCRESULTS(TRACTION,FORCE,PLOT) = ALL
  DISPLACEMENT(PLOT,REAL) = ALL
  EDE(PLOT,AVERAGE,RMAG) = ALL
  FORCE(PLOT,REAL,CENTER) = ALL
  MEFFMASS(PRINT,SUMMARY) = YES
  PRESSURE(PLOT,REAL,TOTAL) = ALL
  SPCFORCES(PLOT,REAL) = ALL
  STRESS(PLOT,REAL,VONMISES,CENTER) = ALL
$*
$* BULK DATA
$*
BEGIN BULK
$*
$* SOLUTION CARDS
$*
TABDMP1      300       G                                                +
+         0.00000.0000001.000+100.020000    ENDT
$ EIGEN METHODS
$*  Modeling Object: Real Eigenvalue - Lanczos
EIGRL, 101, 10.0000, 10.0E3, 50, 0, 7, , MASS
$* Modeling Object: Complex Eigenvalue - Lanczos
EIGC,102, CLAN, , , , , 50
$
$* PARAM CARDS
$*
PARAM      K6ROT100.0000
PARAM     OIBULK     YES
PARAM    OMACHPR     YES
PARAM       POST      -2
PARAM    POSTEXT     YES
PARAM    UNITSYS   MN-MM
$*"""+'\n')

def headers_modes(file):
    #Used for modal analysis
    file.write(r"""
    $*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*                    Modal analysis of a hydrofoil
$*                    for Simcenter Nastran version 2020.2
$*
$*
$ EXECUTIVE CONTROL
$*
ID, modal analysis
SOL 103
CEND
$*
$* CASE CONTROL
$*
ECHO = NONE
SPC = 100
METHOD = 101
OUTPUT
AEROF = ALL
APRESSURE = ALL
DISPLACEMENT(PLOT,REAL) = ALL
EDE(PLOT,AVERAGE,RMAG) = ALL
FORCE(PLOT,REAL,CENTER) = ALL
SPCFORCES(PLOT,REAL) = ALL
STRESS(PLOT,REAL,VONMISES,CENTER) = ALL
$*  Step: Subcase - Static Loads 1
SUBCASE 1
  LABEL = Subcase - Eigenvalue Method 1
  METHOD = 101
  OUTPUT
  AEROF = ALL
  APRESSURE = ALL
  DISPLACEMENT(PLOT,REAL) = ALL
  EDE(PLOT,AVERAGE,RMAG) = ALL
  FORCE(PLOT,REAL,CENTER) = ALL
  SPCFORCES(PLOT,REAL) = ALL
  STRESS(PLOT,REAL,VONMISES,CENTER) = ALL
$*
$* BULK DATA
$*
BEGIN BULK
$*
$* SOLUTION CARDS
$*
$*  Modeling Object: Real Eigenvalue - Lanczos1
EIGRL, 101, 10.0000, 30.0E3, 50, 0, 7, , MASS
$* PARAM CARDS
$*
PARAM      K6ROT100.0000
PARAM     OIBULK     YES
PARAM    OMACHPR     YES
PARAM       POST      -2
PARAM    POSTEXT     YES
PARAM    UNITSYS   MN-MM
$*"""+"\n")

def headers_bending(file):
    #Used for debugging the structural part of this program
    file.write(r"""
    $*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$ Bending analysis of a hydrofoil
$*
$* EXECUTIVE CONTROL
$*
ID, static gravity
SOL 101
CEND
$*
$ CASE CONTROL
$*
ECHO = NONE
SPC = 100
OUTPUT
AEROF = ALL
APRESSURE = ALL
DISPLACEMENT(PLOT,REAL) = ALL
EDE(PLOT,AVERAGE,RMAG) = ALL
FORCE(PLOT,REAL,CENTER) = ALL
SPCFORCES(PLOT,REAL) = ALL
STRESS(PLOT,REAL,VONMISES,CENTER) = ALL
$*  Step: Subcase - Static Loads 1
SUBCASE 1
  LABEL = Subcase - Static Loads 1
  LOAD = 100
$*
$ BULK DATA
$*
BEGIN BULK
$*
$* PARAM CARDS
$*
PARAM      K6ROT100.0000
PARAM     OIBULK     YES
PARAM    OMACHPR     YES
PARAM       POST      -2
PARAM    POSTEXT     YES
PARAM    UNITSYS   MN-MM
$*"""+"\n")

def footers(file):
    #Finishing the bdf file and closing the file
    file.write('$*\n')
    file.write('$*\nENDDATA')
    file.close()