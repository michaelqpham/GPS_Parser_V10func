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
import struct
def binary(num):
    return ''.join(bin(c).replace('0b', '').rjust(8, '0') for c in struct.pack('!f', num))

def dec2hex(n):
    """return the hexadecimal string representation of integer n"""
    return "%X" % n
def hex2dec(s):
    """return the integer value of a hexadecimal string s"""
    return int(s, 16)

'''read data from binary file
fid = fileID, nelements = sizeA, dtype=precision, 
'''

def fread(fid, nelements, dtype):
    if dtype is np.str:
        dt = np.uint8  # WARNING: assuming 8-bit ASCII for np.str!
    else:
        dt = dtype

    data_array = np.fromfile(fid, dt, nelements)
    data_array.shape = (nelements, 1)
    return data_array


def fread1(fid, nelements, dtype, machinefmt):
    dtype.byteorder = machinefmt

    '''
    Order for reading bytes in the file: 
    'ieee-le' = little-endian ordering
    -----
    Specify byte order ">" is little-endian
    machinefmt = "little"
    '''
    data_array = np.fromfile(fid, dtype, nelements)
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
        fid3501 = open(msg3501Filename,'w+')
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

    ofid = open(os.path.join(datadir, summfile), 'w+')
    tfid = open(os.path.join(datadir, parsfile), 'w+')
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
    ofid.write('%s \n' % '======================================================================');
    ofid.write('%s \n\n' % ('   ' + Estampa));
    ofid.write('%s \n' % ('   GPS File: ' + gpsname));
    ofid.write('%s \n' % ('======================================================================'));

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

    while WordTest != hex2dec('81FF'): 
        atByte += 1
        gfid.seek(atByte,0)
        WordTest  = fread(gfid,1, np.uint16)

    GarBytes = np.int64(0)
    count = np.int64(0)
    for line in gfid: # while not eof
        if atByte > count+100000:
            count = atByte
            # add in a waitbar
        tfid.write('%s \n'% '------------------------------')
        tfid.write('%s \n'% ('Byte at start = ' + str(atByte)))
        # absolute positioning
        gfid.seek(atByte+0, 0)
        SynWrd = fread(gfid, 1, np.uint16)[0]
        gfid.seek(atByte+2, 0)
        MssgID = fread(gfid, 1, np.uint16)[0]              # Message 3 ID
        gfid.seek(atByte+4, 0)
        WrdCnt = fread(gfid, 1, np.uint16)[0]              # Word Count
        gfid.seek(atByte+6, 0)
        FlgWrd = fread(gfid, 1, np.uint16)[0]              # Flag Word
        gfid.seek(atByte+8, 0)
        HdrCkSm = fread(gfid, 1, np.uint16)[0]             # Header Check Sum

        # get the bitwise complement
        hdrsumhex = dec2hex(~(np.uint32(SynWrd) + np.uint32(MssgID) + np.uint32(WrdCnt) + np.uint32(FlgWrd)))

        if (len(hdrsumhex) > 4):
            expcthsm = hex2dec(hdrsumhex[(len(hdrsumhex)-4):len(hdrsumhex)])+1
        else:
            expcthsm = hex2dec(hdrsumhex)+1

        if (float(HdrCkSm)-expcthsm) == 0:
            MsgCount += 1
            tfid.write('%s\n'% ['Message Count = ' + str(MsgCount) + ' >> Message Type: ' + str(MssgID) + ' >>> Data Word Count: ' + str(WrdCnt)])
        # TEST FOR A TRUNCATED MESSAGE WHILE IN HEADER OR DATA THEN. IF IT IS, THEN IT WILL
        # NOT PARSE THE MESSAGE DATA.
        # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
        if (atByte + 12 >= sp_bytes):
            TruncEnd = 1
        else:
            TruncEnd = ((atByte + 10 + np.int64(WrdCnt)*2 + 2) > sp_bytes)
        #  RUN DATA CHECKSUM FOR EVERY IDENTIFIED MESSAGE IF THERE IS
        #  CONTENT IN MESSAGE
        # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
        if (WrdCnt != 0) and not TruncEnd:
            DataByte = atByte + 10
            MsgData = np.zeros((1, WrdCnt[0]+1))

            # include np.int64(WrdCnt)
            for iw in range(1, np.int64(WrdCnt)+1):
                gfid.seek(DataByte+2*(iw-1), 0)
                MsgData[0][iw] = fread(gfid, 1, np.uint16)[0]
            gfid.seek(DataByte+2*iw, 0)
            DataChckSum = fread(gfid, 1, np.uint16)[0]
            DataChckSum = float(DataChckSum)
            # convert the decimal to binary and cut out the '0b' at the front
            datsumbin = str(binary(float(np.sum(MsgData))))

            if len(datsumbin) > 16:
                # python's indexing starts at 0 as opposed to matlab's start at 1
                # convert binary to decimal with int(n, 2)
                if ((int((datsumbin[len(datsumbin)-16:len(datsumbin)]), 2)) != 0):
                    expctdsm = (int((datsumbin[len(datsumbin)-16:len(datsumbin)]), 2))
                else:
                    expctdsm = (int((datsumbin[len(datsumbin)-16:len(datsumbin)]), 2)) + (2**16)
            else:
                expctdsm = int(datsumbin, 2)
                if datsumbin != 0:
                    expctdsm += 1
            if (float(DataChckSum) + expctdsm)== 2**16:
                tfid.write('%s\n' % '\tData CheckSum OK')
            else:
                DatCkSmFail += 1
                tfid.write('%s\n' % ' *** ERROR *** Data CheckSum FAILS!')

            # clear the variables. Not strictly needed in python (?)
            del DataByte, MsgData, iw, datsumbin, expctdsm

        # IDENTIFY MESSAGE AND APPROPRIATELY ADD TO ACCOUNTING
        #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
        if MssgID == 3:
            Msg3Count += 1
            tfid.write('%s\n' %('Mssg 3 Count: ' + str(Msg3Count)))
            gfid.seek(atByte + 10, 0)
            # GPStime
            # TODO: Figure out how to specify byteorder
            # python's native byteorder (search with $ sys.byteorder) 
            # is "little" (little-endian)

            GPSd1[:4] = fread(fread, 1, np.dtype('uint16') )[0]
            gtW1bP = bin(GPSd1[1])
            gtW2bP = bin(GPSd1[0])
            gtW3bP = bin(GPSd1[3])
            gtW4bP = bin(GPSd1[2])

    return [MSGCOUNT, REPCOUNT, GPSTimSec]
