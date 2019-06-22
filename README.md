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

## troubleshooting 
### 1. THaSlotData.h:160: int Decoder::THaSlotData::getData(int, int) const: Assertion `chan >= 0 && chan < (int)maxc && hit >= 0 && hit < numHits[chan]' failed. 

This issues is caused by the Analyzer. The Max Chan and the DATA is not large enough to handle the MPD slot more than 16. Edit the file '  hana_decode/THaCrateMap.C'. changed the following line:

    const UShort_t THaCrateMap::MAXCHAN =4096 * 2;   
    const UShort_t THaCrateMap::MAXDATA = 32768 * 2;


## some other issues

### 1. in the mapping the gem_id in the mapping does not match the configuration
