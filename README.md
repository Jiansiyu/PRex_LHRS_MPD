# PREX-counting
Reconstruction code for PREX GEMs in counting mode

This code is is built against the TreeSearch library at $TREESEARCH

https://github.com/JeffersonLab/TreeSearch/

Contains:
    MPDModule
    Decoder for MPD/APV25 used for GEMs

    MPD GEM Classes
    Decoding, clustering, hit-finding based on TreeSearch::GEMPlane

    replay/replay.C
    Example driver script for analysis of GEM data


    DB/db_*
    Example databases for classes
    
## modification 
* change the MPDGEMPlane::Decode. rewrite part of the code , add common mode suppression before Pedestal Suppression. Changed the input value for ChargeDep to samples after CM suppresion and Pedestal Suppression 
* MPDGEMPlane::ChargeDep, disabled the deconverlution to solve the negative adc issues. 
* changed the Pedestals to the UVa standard pedestal 
* some return raw data size does not match '128*Nsamples', ignored those APVs
* change the PREXStand class, add the function used for project the GEM (Transport Coordination System) to the Focal Plane Coordination and to the Target Coordination System.
* Add another file in the database so as to input the Projection Matrix. Latter will need to merge to the main GEM database

## troubleshooting 
#### 1. THaSlotData.h:160: int Decoder::THaSlotData::getData(int, int) const: Assertion `chan >= 0 && chan < (int)maxc && hit >= 0 && hit < numHits[chan]' failed. 

This issues is caused by the Analyzer. The Max Chan and the DATA is not large enough to handle the MPD slot more than 16. Edit the file '  hana_decode/THaCrateMap.C'. changed the following line:

    const UShort_t THaCrateMap::MAXCHAN =4096 * 2;   
    const UShort_t THaCrateMap::MAXDATA = 32768 * 2;


## some other issues and TODOs

* 1. in the mapping the gem_id in the mapping does not match the configuration
* 2. The formula for the Projection need to check again. Maybe not very right
* 3. Add the rotation correction algorithm and also new Eular angle correction parameters need to be added in the database 



## branch instructions

##### 1. prex_rhrs

Main branch for the PRex_MPD. The mainly test branch for the development.

##### 2. prex_rhrs_stable

Main branch for the PRex_MPD. It should be the most updated stable version. 

