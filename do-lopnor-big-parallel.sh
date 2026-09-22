#!/bin/bash
##
##  "DO"-script for a Radiative3D run:
##
##  Copy and edit this file to easily specify (and remember)
##  parameters used for a particular run.
##
##  If $1 is set, it is used as an extra identifier token in the
##  output directory name.
##
##  (Note: Do not edit line 3.  It is used by subsequent scripts to
##  identify this file as a "do script". (See wikify.sh.))
##
source scripts/do-fundamentals.sh
##
##  This DO-SCRIPT exercises shared-memory parallel simulation with the
##  LOP NOR cylinder model while retaining the model's isotropic scattering.
##

## One-liner description: (Keep this BRIEF.)
##
INTENT="Reproducible four-worker Lop Nor waveform simulation."
CAMPAIGN="Parallel capability example"

SIMTARGET="waveform"          # Choice: 'waveform' or 'video'. Affects
                              # defaults not otherwise specified.
MODIDX=1  # LOP NOR           # Model Index: Selects from custom coded models.
                              # 1: Lop Nor (base or moho depends on COMPARGS),
                              # 5: North Sea Crust Pinch model,
                              # 8: Crust Upthrust model
                              # 40: Halfspace model

event=eq        # 'eq' or 'expl' - See case statement below for choices
EQISOFRAC=0.0   # Iso fraction for EQ event (choose from range [-1.0, 1.0])
                # (Or choose as Iso Angle in range [-90, -1), (1, 90].)

##
##  Lop Nor geography and event source parameters:
##

XINXYZ=425.54,-169.53,0.98    # Lop Nor, surface, site of Xinjiang Quake 030313
MAKXYZ=-102.27,430.84,0.60   # Station MAK
WUSXYZ=-390.04,-167.18,1.457 # Station WUS

FREQ=2.0                      # Phonon frequency to model
NUMPHONS=10M                  # Same workload as do-lopnor-big.sh
RECTIME=600                   # Recording duration of seismometers (seconds).
BINSIZE=2.00                  # Seismometer time-bin size in seconds
GATHER=40.0                   # Terminal gather radius, in kilometers.
CYLRAD=1200                   # Total (cylindrical) radius of model.
                              # (Phonons exceeding this radius from XY
                              # = (0,0), or this travel time from t=0,
                              # are abandoned.)

SCAT1=0.8,0.01,0.5,0.2,50     # Scat Args (nu,eps,a,kappa) and Q in sedi's
SCAT2=0.8,0.01,0.5,0.3,1000   # Scat Args (nu,eps,a,kappa) and Q in crust
SCAT3=0.8,0.01,0.7,0.5,300    # Scat Args (nu,eps,a,kappa) and Q in mantle
                              # (Q values specified are Q_s values.)

COMPARGS=$SCAT1,$SCAT2,$SCAT3

FLATTEN="--flatten"          # Apply Earth-flattening transformation.

# Explicit worker count and seed make this a reproducible parallel recipe.
# No scattering-length override is supplied: the compiled model remains
# isotropic, providing a clean parallelism example.
ADDITIONAL="--workers=4 --seed=0x5eedc0de12345678"

case "$event" in
    expl)   # Generic explosion
        SOURCELOC=425.54,-169.53,-1.02
        SOURCETYP=EXPL
        ;;
    eq)     # Xinjiang earthquake
        SOURCELOC=425.54,-169.53,-31.02
        SOURCETYP=SDR,125,40,90,$EQISOFRAC
        ;;
    *)
        echo Event code not recognized.
        exit
        ;;
esac

SEISORIG1=$XINXYZ
SEISORIG2=$XINXYZ
SEISDEST1=$WUSXYZ
SEISDEST2=$MAKXYZ
SEIS1=--seis-p2p=$SEISORIG1,$SEISDEST1,1.0,2.0,$GATHER,160
SEIS2=--seis-p2p=$SEISORIG2,$SEISDEST2,1.0,2.0,$GATHER,160

## RUN THE SIMULATION:
##

PopDefaults $SIMTARGET        ## Defined in do-fundamentals.sh
CheckBuildStatus              ##  ''
CreateOutputDirectory $@      ##  ''
RunSimulation || exit $?      ##  ''
