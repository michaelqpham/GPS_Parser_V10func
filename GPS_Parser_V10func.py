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
import math
import datetime as dt 
import numpy as np 

def dec2hex(n):
    """return the hexadecimal string representation of integer n"""
    return "%X" % n
def hex2dec(s):
    """return the integer value of a hexadecimal string s"""
    return int(s, 16)

'''read data from binary file'''
def fread(fid, nelements, dtype):
    if dtype is np.str:
     dt = np.uint8  # WARNING: assuming 8-bit ASCII for np.str!
    else:
     dt = dtype

    data_array = np.fromfile(fid, dt, nelements)
    data_array.shape = (nelements, 1)

    return data_array

def GPS_Parser_V10func(filename, datadir, GPSFigs, PAR3501):
    Vnum = 'V10'
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
    gfid = open(gpsname, 'rb')

    '''
    %-------------------------------------------------------------------------------------
    % DETERMINE IF MSG3501 ALREADY EXISTS AND IF NOT, CREATE IT.
    %-------------------------------------------------------------------------------------
    '''
    msg3501Filename = gpsname.replace('.bin','') + '-Msg3501.csv'
    parsed3501 = os.path.exists(msg3501Filename)

    if PAR3501 and not parsed3501:
        fid3501 = open(msg3501Filename,'wb')
        fid3501.write('%s \n' % ('%============================================================================='))
        # TODO Brackets were around (GPS_Parser_V10func.m, line 101)
        fid3501.write('%s \n' % ('%NGDCS Check ' + Vnum + ': Message 3501 Navigation Solution of ' + gpsname))
        fid3501.write('%s \n' % ('%* If present "Date" and "Abs Time" are derived and not original Msg3501 data'))
        fid3501.write('%s \n' % ('%============================================================================='))
        fid3501.write('%s \n' % ('%Msg3501#, StartByte, Flag, Time, Lat, Lon, Alt, VelN, VelE, VelUp, Pitch, Roll, HeadTrue, DataCheckSum'))
    '''
    %-------------------------------------------------------------------------------------
    % READ C-MIGITS III DATA AND FIND MESSAGE ID3 FOR EACH PULSE. BUT FIRST,
    % SEEK THE FIRST hx81FF CODE FOR THE FIRST COMPLETE MESSAGE.
    %-------------------------------------------------------------------------------------
    '''

    allfile  = filename + '-parsed.txt'
    summfile = filename + '-summtemp.txt'
    parsfile = filename + '-parsetemp.txt'

    ofid = open(os.path.join(datadir, summfile), 'wb')
    tfid = open(os.path.join(datadir, parsfile), 'wb')
    sp_bytes = np.int64(os.path.getsize(os.path.join(datadir,gpsname)))
    nwrd = np.int64(math.ceil(sp_bytes / 2.0))
    nams = np.int64(nwrd/5.0)

    # print("sp: ", sp_bytes)
    # print("nwrd: ", nwrd)
    # print("nams: ", nams)

    GPSTimSec = np.zeros(nams)
    GPSTimSec = np.zeros(nams)
    MSGPCntrs = np.zeros(5, nams)

    Estampa = 'GPS_Parser_' + Vnum +' on ' + str(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    # print(Estampa)
    ofid.write('%s \n' % ('======================================================================'))
    ofid.write('%s \n\n' % ('   ' + Estampa))
    ofid.write('%s \n' % ('   GPS File: ' + gpsname))
    ofid.write('%s \n' % ('======================================================================'))

    atByte = np.int64(0)
    MsgCount = np.int64(0)
    Msg3Count = np.int64(0)
    Msg3500Count = np.int64(0)
    Msg3501Count = np.int64(0)
    Msg3502Count = np.int64(0)
    Msg3512Count = np.int64(0)
    Msg3623Count = np.int64(0)
    Msg3500PartCount = np.int64(0)
    Msg3501PartCount = np.int64(0)
    Msg3502PartCount = np.int64(0)
    Msg3512PartCount = np.int64(0)
    Msg3623PartCount = np.int64(0)
    Msg3500CurrMode  = np.int64(0)
    pidx = 1
    DatCkSmFail = 0    
    gfid.seek(0, 0) # start at beginning of file
    WordTest = fread(gfid, 1, np.uint16) # read the file into 16-bit unsigned ints
    
    
    return [MSGCOUNT, REPCOUNT, GPSTimSec]
