'''
%[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]
%
% NGDCS CHECK-CODE SUITE FOR IMAGING SPECTROMETERS - NTR#49120
%
% Copyright 2016, by the California Institute of Technology. ALL
% RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
%
% This software may be subject to U.S. export control laws and regulations.
% By accepting this document, the user agrees to comply with all applicable
% U.S. export laws and regulations.  User has the responsibility to obtain
% export licenses, or other export authority as may be required before
% exporting such information to foreign countries or providing access to
% foreign persons.
%
% Version 10
% Daniel Nunes, Didier Keymeulen, Jacky Lee, JPL/Caltech, Michaek Pham
% DCN, 2015/FEB/26
% Daniel.Nunes@jpl.nasa.gov
%
%--------------------------------------------------------------------------
% Matlab function designed to read CMIGITSII or SND500 navigation data.
%
% INPUT
%   1 - filename : name of binary file containing GPS data. It can have
%                suffixes "gps", "nav" or "cmigits"
%
%   2 - datadir  : path to the filename
%
%   3 - GPSFigs  : true/false flag, which when set to true will create 
%                figures of the GPS analysis
%
%   4 - PAR3501  : true/false flag, which when set to true will generate a
%                *.csv file containing the data for all Msg3501 messages
%
% OUTPUT
%   1 - MSGCOUTN : structure with count of Messages (0003,3500,3501,3502,
%                  3512,3623) identified in GPS file
%
%   2 - REPCOUNT : structure with counts of total bytes, unassigned bytes,
%                  and failed checksums
%
%   3 - GPSTimSec: array containing the detected timestamps (in sec) 
%                  furnished bythe GPS/INS unit in the GPS file
%
% EXTERNAL SUBROUTINES & FUNCTIONS
%   1) generateAbsTimeFunc - takes a user specified Msg3501.csv and 
%                            generate a new file named Msg3501AbsTime.csv
%           
%                   USAGE: generateAbsTimeFunc( datadir, msg3501Filename )
%
%   1.1)toAbsTimeFunc - returns an absolute time of the week after 
%                        converting the relative time in decimal seconds 
%                        from a the 3501 message
%
%                   USAGE:   dayAndTime = toAbsTimeFunc( timeRelative )
%
%
%[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]
'''
import os.path # Import os.path to search file paths

def GPS_Parser_V10func(filename, datadir, GPSFigs, PAR3501):
    '''
    %=====================================================================================
    % U S E R    P A R A M E T E R S
    %=====================================================================================
    '''
    Vnum = 'V10'
    PssblMsg = [3, 3500, 3501, 3502, 3512, 3623]
    GPSModes = ['(1) Test ','(2) Initialization ','(3) *not used* ',
                '(4) Fine Alignment ','(5) Air Alignment ','(6) Transfer Alignment ',
                '(7) Air Navigation ','(8) *not used* ','(9) GPS Only ']
    # MSGCOUNT = struct('M0003',[],'M3500',[],'M3501',[],'M3502',[],'M3512',[],'M3623',[])
    # REPCOUNT = struct('TotalBytes',[],'UnassignedBytes',[],'FailedChkSum',[],
    #                   'NoGPSModes',[],'GPSModes',[],'PctNavMode',[])
    MSGCOUNT = {'M0003':[], 'M3500':[], 'M3501':[], 'M3502':[], 'M3512':[], 'M3623':[]}
    REPCOUNT = {'TotalBytes': [], 'UnassignedBytes':[], 'FailedChkSum':[], 'NoGPSModes':[], 
        'GPSModes':[], 'FailedChkSum':[]}
    '''
    %=====================================================================================

    %-------------------------------------------------------------------------------------
    % INITIAL USER INTERACTIONS - CHOOSE FILE AND OPT FOR FIGURE
    %-------------------------------------------------------------------------------------
    '''
    # hgps = waitbar(0,['Processing GPS data...'])

    '''
    %-------------------------------------------------------------------------------------
    % OPEN THE GPS FILE
    %-------------------------------------------------------------------------------------
    '''
    gpsname = filename
    gfid = open(gpsname, 'r')

    '''
    %-------------------------------------------------------------------------------------
    % DETERMINE IF MSG3501 ALREADY EXISTS AND IF NOT, CREATE IT.
    %-------------------------------------------------------------------------------------
    '''
    msg3501Filename = gpsname.replace('.bin','') + '-Msg3501.csv'
    parsed3501 = os.path.exists(msg3501Filename)

    if PAR3501 and not parsed3501:
        fid3501 = open(msg3501Filename,'w')
        fid3501.write('%s \n', '%=============================================================================');
        # TODO Brackets were around (GPS_Parser_V10func.m, line 101)
        fid3501.write('%s \n', '%NGDCS Check ' + Vnum + ': Message 3501 Navigation Solution of ' + gpsname);
        fid3501.write('%s \n', '%* If present "Date" and "Abs Time" are derived and not original Msg3501 data');
        fid3501.write('%s \n', '%=============================================================================');
        fid3501.write('%s \n', 'Msg3501#, StartByte, Flag, Time, Lat, Lon, Alt, VelN, VelE, VelUp, Pitch, Roll, HeadTrue, DataCheckSum');

    return [MSGCOUNT, REPCOUNT, GPSTimSec]
