#Generating the headers for a BDF file
#Author: Danick Lamoureux
#LM2 project under Frédérick Gosselin's supervision
#Date: 2022-05-05

def headers_hydro(file, P1, nmodes = 10):
    #Headers to use for full hydroelastic analysis. Requires modal and coupled acoustics analysis first


    file.write(r"""$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*                    Aeroelasticity approximation of hydrodynamic damping
$*                    for Simcenter Nastran version 2020.2
$*
$*        ANALYSIS TYPE: Structural
$*        SOLUTION NAME: Flutter
$*        SOLUTION TYPE: SOL 145 Aeroelastic Flutter
$*
$*                UNITS: mm (milli-newton)
$*                      ... LENGTH : mm
$*                      ... TIME   : sec
$*                      ... MASS   : kilogram (kg)
$*                      ... TEMPERATURE : deg Celsius
$*                      ... FORCE  : milli-newton
$*                      ... THERMAL ENERGY : mN-mm (micro-joule)
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* FILE MANAGEMENT
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* EXECUTIVE CONTROL
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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
    #À garder pour le calcul en cascades
    # for i in range(nmodes):
    #     file.write("\nMATGEN ,/QQQ"+str(i+1)+"/1/1 $\n")
    #     file.write("ADD QQQ"+str(i+1)+",/QQ"+str(i+1)+"/"+str(P2[i,i])+" $")
    #     if i+1>1:
    #         if i+1==2:
    #             file.write("\nMERGE QQ1,ZEROT1,ZERO1,QQ2,CP1,/B1 $")
    #         else:
    #             file.write("\nMERGE B"+str(i-1)+",ZEROT"+str(i)+",ZERO"+str(i)+",QQ"+str(i+1)+",CP"+str(i)+",/B"+str(i)+" $")
    # for i in range(nmodes):
    #     file.write("\nMATGEN ,/RRR"+str(i+1)+"/1/1 $\n")
    #     file.write("ADD RRR"+str(i+1)+",/RR"+str(i+1)+"/"+str(P3[i,i])+" $")
    #     if i+1>1:
    #         if i+1==2:
    #             file.write("\nMERGE RR1,ZEROT1,ZERO1,RR2,CP1,/C1 $")
    #         else:
    #             file.write("\nMERGE C"+str(i-1)+",ZEROT"+str(i)+",ZERO"+str(i)+",RR"+str(i+1)+",CP"+str(i)+",/C"+str(i)+" $")
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
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* CASE CONTROL
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* BULK DATA
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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

def headers_simulated_modes(file):
    #Headers to use for full hydroelastic analysis. Requires modal and coupled acoustics analysis first
    file.write(r"""$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*                    Aeroelasticity approximation of hydrodynamic damping
$*                    for Simcenter Nastran version 2020.2
$*
$*        ANALYSIS TYPE: Structural
$*        SOLUTION NAME: Flutter
$*        SOLUTION TYPE: SOL 145 Aeroelastic Flutter
$*
$*                UNITS: mm (milli-newton)
$*                      ... LENGTH : mm
$*                      ... TIME   : sec
$*                      ... MASS   : kilogram (kg)
$*                      ... TEMPERATURE : deg Celsius
$*                      ... FORCE  : milli-newton
$*                      ... THERMAL ENERGY : mN-mm (micro-joule)
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* FILE MANAGEMENT
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* EXECUTIVE CONTROL
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
ID, Modal analysis, added mass
SOL 107
CEND
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* CASE CONTROL
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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
  M2GG = MAJ
  K2GG = KAJ
  B2GG = BAJ
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* BULK DATA
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
BEGIN BULK
$*
$* SOLUTION CARDS
$*
$*  Modeling Object: Real Eigenvalue - Lanczos1
EIGRL, 101, 10.0000, 10.0E3, , 0, 7, , MASS
$* Modeling Object: Complex Eigenvalue - Lanczos
EIGC,102, CLAN, , , , , 20
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
$*                    Aeroelasticity approximation of hydrodynamic damping
$*                    for Simcenter Nastran version 2020.2
$*
$*        ANALYSIS TYPE: Structural
$*        SOLUTION NAME: Flutter
$*        SOLUTION TYPE: SOL 145 Aeroelastic Flutter
$*
$*                UNITS: mm (milli-newton)
$*                      ... LENGTH : mm
$*                      ... TIME   : sec
$*                      ... MASS   : kilogram (kg)
$*                      ... TEMPERATURE : deg Celsius
$*                      ... FORCE  : milli-newton
$*                      ... THERMAL ENERGY : mN-mm (micro-joule)
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* FILE MANAGEMENT
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* EXECUTIVE CONTROL
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
ID, Aeroelastic analysis, no added mass
SOL 145
CEND
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* CASE CONTROL
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* BULK DATA
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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
$*                    Aeroelasticity approximation of hydrodynamic damping
$*                    for Simcenter Nastran version 2020.2
$*
$*
$*
$*        ANALYSIS TYPE: Vibro-Acoustic
$*        SOLUTION NAME: Solution 1
$*        SOLUTION TYPE: SOL 103 Real Eigenvalues
$*
$*
$*                UNITS: mm (milli-newton)
$*                      ... LENGTH : mm
$*                      ... TIME   : sec
$*                      ... MASS   : kilogram (kg)
$*                      ... TEMPERATURE : deg Celsius
$*                      ... FORCE  : milli-newton
$*                      ... THERMAL ENERGY : mN-mm (micro-joule)
$*
$* IMPORTANT NOTE:
$*     This banner was generated by Simcenter and altering this 
$*     information may compromise the pre and post processing of results
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* FILE MANAGEMENT
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* EXECUTIVE CONTROL
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
ID, Acoustics analysis, no added damping, uncoupled fluid-structure
SOL 103
CEND
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* CASE CONTROL
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* BULK DATA
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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
$*                    Aeroelasticity approximation of hydrodynamic damping
$*                    for Simcenter Nastran version 2020.2
$*
$*
$*
$*        ANALYSIS TYPE: Vibro-Acoustic
$*        SOLUTION NAME: Solution 1
$*        SOLUTION TYPE: SOL 107 Complex Eigenvalues
$*
$*
$*                UNITS: mm (milli-newton)
$*                      ... LENGTH : mm
$*                      ... TIME   : sec
$*                      ... MASS   : kilogram (kg)
$*                      ... TEMPERATURE : deg Celsius
$*                      ... FORCE  : milli-newton
$*                      ... THERMAL ENERGY : mN-mm (micro-joule)
$*
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* FILE MANAGEMENT
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* EXECUTIVE CONTROL
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
ID, Acoustics analysis, no added damping, coupled fluid-structure
SOL 107
CEND
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* CASE CONTROL
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* BULK DATA
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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

def headers_acoustics_complex_EXTOUT(file, nsol, nflu):
    #Coupled acoustics that print out the MAJ in OP2 format
    file.write(r"""$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*                    Aeroelasticity approximation of hydrodynamic damping
$*                    for Simcenter Nastran version 2020.2
$*
$*
$*
$*        ANALYSIS TYPE: Vibro-Acoustic
$*        SOLUTION NAME: Solution 1
$*        SOLUTION TYPE: SOL 103 Real Eigenvalues
$*
$*
$*                UNITS: mm (milli-newton)
$*                      ... LENGTH : mm
$*                      ... TIME   : sec
$*                      ... MASS   : kilogram (kg)
$*                      ... TEMPERATURE : deg Celsius
$*                      ... FORCE  : milli-newton
$*                      ... THERMAL ENERGY : mN-mm (micro-joule)
$*
$* IMPORTANT NOTE:
$*     This banner was generated by Simcenter and altering this 
$*     information may compromise the pre and post processing of results
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* FILE MANAGEMENT
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* EXECUTIVE CONTROL
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
ID, Acoustics analysis, no added damping, coupled fluid-structure
SOL 107
$COMPILE PHASE1DR NOLIST NOREF $
$MALTER '$MALTER:(KGG, BGG, MGG, K4GG, PG)' $
$MALTER '$MALTER:STATIC AND DYNAMIC (KAA, MAA, BAA, K4AA, PA)' $
$MATPRN BGG $
$COMPILE PHASE1DR NOLIST NOREF $
$MALTER '$MALTER:STATIC AND DYNAMIC (KAA, MAA' $
$MATGEN /CP/6/"""+str(nsol+nflu)+'/'+str(nsol)+'/'+str(nflu)+""" $
$MATGEN /RP/6/"""+str(nsol+nflu)+'/0/'+str(nsol)+'/'+str(nflu)+""" $
$PARTN AGG,CP,RP/,,,ACSF/1/0/ $
$PARTN KGG,CP,RP/,,KF,/1/0/ $
$TRNSP ACSF/ACSFT $
$SOLVE KF,ACSFT,,,/X $
$MPYAD ACSF,X,/MAJ/////6$
$OUTPUT2  MAJ,,,,//0/71///'MATRIX'$
$MATPRN MAJ// $
CEND
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* CASE CONTROL
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
ECHO = NONE
PARAM, EXTOUT, DMIGPCH
$SPC = 100
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
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* BULK DATA
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
BEGIN BULK
$*
$* SOLUTION CARDS
$*
TABDMP1      300       G                                                +
+         0.00000.0000001.000+100.020000    ENDT
$ EIGEN METHODS
$*  Modeling Object: Real Eigenvalue - Lanczos
EIGRL, 101, 10.0000, 10.0E3, 20, 0, 7, , MASS
$* Modeling Object: Complex Eigenvalue - Lanczos
EIGC,102, CLAN, , , , , 20
$
$* PARAM CARDS
$*
PARAM      K6ROT100.0000
PARAM     OIBULK     YES
PARAM    OMACHPR     YES
PARAM       POST      -2
PARAM    POSTEXT     YES
PARAM    UNITSYS   MN-MM
PARAM, EXTOUT, DMIGPCH
$*"""+'\n')

def headers_modes(file):
    #Used for modal analysis
    file.write(r"""
    $*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*                    Aeroelasticity approximation of hydrodynamic damping
$*                    for Simcenter Nastran version 2020.2
$*
$*        ANALYSIS TYPE: Structural
$*        SOLUTION NAME: Modes
$*        SOLUTION TYPE: SOL 103 Real Eigenvalues
$*
$*                UNITS: mm (milli-newton)
$*                      ... LENGTH : mm
$*                      ... TIME   : sec
$*                      ... MASS   : kilogram (kg)
$*                      ... TEMPERATURE : deg Celsius
$*                      ... FORCE  : milli-newton
$*                      ... THERMAL ENERGY : mN-mm (micro-joule)
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* FILE MANAGEMENT
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* EXECUTIVE CONTROL
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
ID, modal analysis
SOL 103
CEND
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* CASE CONTROL
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* BULK DATA
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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

def headers_modes_EXTOUT(file):
    #Used to print the MAJ in DMIG cards
    file.write(r"""
    $*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*                    Aeroelasticity approximation of hydrodynamic damping
$*                    for Simcenter Nastran version 2020.2
$*
$*        ANALYSIS TYPE: Structural
$*        SOLUTION NAME: Modes
$*        SOLUTION TYPE: SOL 103 Real Eigenvalues
$*
$*                UNITS: mm (milli-newton)
$*                      ... LENGTH : mm
$*                      ... TIME   : sec
$*                      ... MASS   : kilogram (kg)
$*                      ... TEMPERATURE : deg Celsius
$*                      ... FORCE  : milli-newton
$*                      ... THERMAL ENERGY : mN-mm (micro-joule)
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* FILE MANAGEMENT
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* EXECUTIVE CONTROL
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
ID, modal analysis
SOL 103
CEND
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* CASE CONTROL
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
ECHO = NONE
PARAM, EXTOUT, DMIGPCH
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
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* BULK DATA
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
BEGIN BULK
$*
$* SOLUTION CARDS
$*
$*  Modeling Object: Real Eigenvalue - Lanczos1
EIGRL, 101, 10.0000, 10.0E3, 20, 0, 7, , MASS
$* PARAM CARDS
$*
PARAM      K6ROT100.0000
PARAM     OIBULK     YES
PARAM    OMACHPR     YES
PARAM       POST      -2
PARAM    POSTEXT     YES
PARAM    UNITSYS   MN-MM
$*"""+"\n")

def headers_addedmass_modes_EXTOUT(file):
    #Used to print the MAJ in DMIG cards
    file.write(r"""
    $*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*                    Aeroelasticity approximation of hydrodynamic damping
$*                    for Simcenter Nastran version 2020.2
$*
$*        ANALYSIS TYPE: Structural
$*        SOLUTION NAME: Modes
$*        SOLUTION TYPE: SOL 103 Real Eigenvalues
$*
$*                UNITS: mm (milli-newton)
$*                      ... LENGTH : mm
$*                      ... TIME   : sec
$*                      ... MASS   : kilogram (kg)
$*                      ... TEMPERATURE : deg Celsius
$*                      ... FORCE  : milli-newton
$*                      ... THERMAL ENERGY : mN-mm (micro-joule)
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* FILE MANAGEMENT
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* EXECUTIVE CONTROL
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
ID, modal analysis
SOL 103
CEND
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* CASE CONTROL
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
ECHO = NONE
PARAM, EXTOUT, DMIGPCH
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
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* BULK DATA
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
BEGIN BULK
$*
$* SOLUTION CARDS
$*
$*  Modeling Object: Real Eigenvalue - Lanczos1
EIGRL, 101, 10.0000, 10.0E3, 20, 0, 7, , MASS
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
$*                    Aeroelasticity approximation of hydrodynamic damping
$*                    for Simcenter Nastran version 2020.2
$*
$*        ANALYSIS TYPE: Structural
$*        SOLUTION NAME: Acoustics
$*        SOLUTION TYPE: SOL 101 Static bending
$*
$*                UNITS: mm (milli-newton)
$*                      ... LENGTH : mm
$*                      ... TIME   : sec
$*                      ... MASS   : kilogram (kg)
$*                      ... TEMPERATURE : deg Celsius
$*                      ... FORCE  : milli-newton
$*                      ... THERMAL ENERGY : mN-mm (micro-joule)
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* FILE MANAGEMENT
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* EXECUTIVE CONTROL
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
ID, static gravity
SOL 101
CEND
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* CASE CONTROL
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$* BULK DATA
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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