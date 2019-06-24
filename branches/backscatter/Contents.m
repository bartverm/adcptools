% ADCPtools is a set of functions and classes to process Acoustic Doppler
% Current profiler data
%
%% Reading data:
%   readADCP         - Reads an RDI PD0 binary data file
%   readDeployment   - Read all files in a Winriver deployment
%   readViseaExtern  - Reads Visea extern files
%   readTfiles       - Reads old Transect files
%   readNMEA         - Reads NMEA data from text files
%   readNMEAADCP     - Reads and matches NMEA files with adcp data
%
%% Low-level data processing
%   corADCP          - Transform coordinate system of ADCP velocity data
%   filterADCP       - Filter velocity data and transform to double (m/s)
%   utmADCP          - Obtain ADCP positions in UTM projected coord system
%   depthADCP        - Position of bed detection in each beam wrt ADCP
%   mapADCP          - Position of velocity measurements wrt ADCP 
%   getADCPHeading   - Get ADCP heading, internally or externally measured
%   getExtMisalign   - Misalignment of external heading device and beam 3
%   getGPSvel        - Vessel velocity from GPS positions
%   svADCP           - Computes backscatter strength (dB)
%
%% High-level data processing
%   compQ            - Estimate discharge
%   procTrans        - Processes repeat transect vessel mounted ADCP data
%   tayl_predict     - Defines taylor expansions for use with procTrans
%   tayl_expand      - Defines taylor expansions for use with procTrans
%
%% Packages
%   acoustics        - acoustics package
%   helpers          - helper functions and classes
%
%% Helper functions (usually not called directly)
%   defineNMEA       - This file defines all regular expressions for nmea strings
%   geo2utm          - Transform WGS84 to UTM coordinate system
%   isADCPstruct     - checks whether a variable is an ADCP structure
%   matmult          - Matrix multiplication
%   nmeachecksum     - Compute NMEA checksum
%   readDBS          - readDBT(DBT) interprets a NMEA DBT string
%   readDBT          - readDBT(dbt) interprets a NMEA DBT string
%   readGGA          - readGGA(gga) interprets a NMEA GGA string
%   readGBS          - readGBS(gbs) interprets the csi propietary NMEA message GBS
%   readGLL          - readGLL(gll) interprets a NMEA GLL string
%   readGSA          - readGSA(gsa) interprets a NMEA GSA string
%   readHDT          - readHDT(hdt) interprets a NMEA HDT string
%   readHPR          - readHPR(hpr) interprets a CSI propietary NMEA string $PSAT,HPR
%   readRDENS        - readRDENS(gga) interprets the RDI propietary NMEA message RDENS
%   readRMC          - readRMC(rmc) interprets a NMEA RMC string
%   readVTG          - readVTG(vtg) interprets a NMEA VTG string
%   readZDA          - readZDA(zda) interprets a NMEA ZDA string
%   textscan_checked - To debug textscan errors
