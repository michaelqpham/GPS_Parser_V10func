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
import numpy
import datetime
import numpy.matlib

def GPS_Parser_V10func(filename, datadir, GPSFigs, PAR3501):
    return [MSGCOUNT, REPCOUNT, GPSTimSec]
'''
%=====================================================================================
% U S E R    P A R A M E T E R S
%=====================================================================================
'''
Vnum = 'V10'
PssblMsg = [3, 3500, 3501, 3502, 3512, 3623]
GPSModes = ['(1) Test ','(2) Initialization ','(3) *not used* ',\
            '(4) Fine Alignment ','(5) Air Alignment ','(6) Transfer Alignment ',\
            '(7) Air Navigation ','(8) *not used* ','(9) GPS Only ']
MSGCOUNT = struct('M0003',[],'M3500',[],'M3501',[],'M3502',[],'M3512',[],'M3623',[])
REPCOUNT = struct('TotalBytes',[],'UnassignedBytes',[],'FailedChkSum',[],\
                  'NoGPSModes',[],'GPSModes',[],'PctNavMode',[])
'''
%=====================================================================================

%-------------------------------------------------------------------------------------
% INITIAL USER INTERACTIONS - CHOOSE FILE AND OPT FOR FIGURE
%-------------------------------------------------------------------------------------
'''
# hgps = waitbar(0,['Processing GPS data...'])
#I HAVE NO IDEA HOW TO DO THIS

'''
%-------------------------------------------------------------------------------------
% OPEN THE GPS FILE
%-------------------------------------------------------------------------------------
'''
gpsname = filename
gfid = open(os.path.join(datadir, gpsname), 'r')

'''
%-------------------------------------------------------------------------------------
% DETERMINE IF MSG3501 ALREADY EXISTS AND IF NOT, CREATE IT.
%-------------------------------------------------------------------------------------
'''
msg3501Filename = gpsname.replace('.bin','') + '-Msg3501.csv'
parsed3501 = os.path.exists(msg3501Filename)

if PAR3501 and not parsed3501:
    fid3501 = open(os.path.join(datadir, msg3501Filename),'w')
    fid3501.write('%s \n', '%=============================================================================')
    fid3501.write('%s \n',['%NGDCS Check ' + Vnum + ': Message 3501 Navigation Solution of ' + gpsname])
    fid3501.write('%s \n', '%* If present "Date" and "Abs Time" are derived and not original Msg3501 data')
    fid3501.write('%s \n', '%=============================================================================')
    fid3501.write('%s \n', '%Msg3501#, StartByte, Flag, Time, Lat, Lon, Alt, VelN, VelE, VelUp, Pitch, Roll, HeadTrue, DataCheckSum')

'''
%-------------------------------------------------------------------------------------
% READ C-MIGITS III DATA AND FIND MESSAGE ID3 FOR EACH PULSE. BUT FIRST,
% SEEK THE FIRST hx81FF CODE FOR THE FIRST COMPLETE MESSAGE.
%-------------------------------------------------------------------------------------
'''

allfile  = filename + '-parsed.txt'
summfile = filename + '-summtemp.txt'
parsfile = filename + '-parsetemp.txt'

ofid = open(os.path.join(datadir,summfile), 'w')
tfid = open(os.path.join(datadir,parsfile), 'w')

sp = os.listdir(os.path.join(datadir,gpsname))
nwrd = sp.nbytes / 2.0
nams = (nwrd/5.0), dtype = numpy.int64
# 
#
#
#I have no idea how to convert to 64 bit int
#

GPSTimSec = numpy.zeros(nams)
GPSTimSec3512 = numpy.zeros(nams)
MsgPCntrs = numpy.zeros(5,nams)

Estampa = ['GPS_Parser_' + Vnum + '  on ' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")]
ofid.write('%s \n','======================================================================')
ofid.write('%s \n\n', ['   ' + Estampa])
ofid.write('%s \n',['   GPS File: ' + gpsname])
ofid.write('%s \n','======================================================================')

'''
NO IDEA HOW TO CHANGE INTO 64 BIT INTEGER

'''
atByte = int64(0)
MsgCount = int64(0)
Msg3Count = int64(0)
Msg3500Count = int64(0)
Msg3501Count = int64(0)
Msg3502Count = int64(0)
Msg3512Count = int64(0)
Msg3623Count = int64(0)
Msg3500PartCount = int64(0)
Msg3501PartCount = int64(0)
Msg3502PartCount = int64(0)
Msg3512PartCount = int64(0)
Msg3623PartCount = int64(0)
Msg3500CurrMode  = int64(0)
pidx = 1
DatCkSmFail = 0
'''

NO IDEA HOW TO CHANGE INTO 64 BIT INTEGER

'''
def fread(fid, nelements, dtype):
     if dtype is numpy.str:
         dt = np.uint8  # WARNING: assuming 8-bit ASCII for np.str!
     else:
         dt = dtype

     data_array = numpy.fromfile(fid, dt, nelements)
     data_array.shape = (nelements, 1)

     return data_array
gfid.seek(0,0) #start at beginning of file
WordTest = fread(gfid, 1, *uint16)
#NOT SURE IF hex2dec is python-ok
while WordTest != hex2dec('81FF'): 
    atByte += 1
    gfid.seek(atByte,0)
    WordTest  = fread(gfid,1, *uint16)

GarBytes = int64(0) #don't know if this works
count = int64(0) #don't know if this works
line = gfid.readline()
while line != '':
    if atByte > count+100000:
        count = atByte
        #%pause(1)
        #waitbar(double(count)/sp.bytes,hgps,['Processing GPS data...' num2str(100*count/sp.bytes) '%']);
        #I HAVE NO IDEA HOW TO DO THIS
    tfid.write('%s \n','------------------------------')
    tfid.write('%s \n',['Byte at start = ' + str(atByte)])
    line = obj.readline()
    
    gfid.seek(atByte+0, 0)
    SynWrd  = fread(gfid,1, *uint16)
    gfid.seek(atByte+2,0)
    MssgID  = fread(gfid,1, *uint16)             # Message 3 ID
    gfid.seek(atByte+4,0)
    WrdCnt  = fread(gfid,1, *uint16)             # Word count
    fseek(gfid,atByte+6, 0)  
    FlgWrd  = fread(gfid,1, *uint16)             # Flag word
    gfid.seek(atByte+8,0)  
    HdrCkSm = fread(gfid,1, *uint16)             # Header Check Sum
    
    #I HAVE NO IDEA HOW TO DO THAT
    hdrsumhex = dec2hex(bitcmp(uint32(SynWrd)+uint32(MssgID)+uint32(WrdCnt)+uint32(FlgWrd)))
    #I HAVE NO IDEA HOW TO DO THAT
    if len(hdrsumhex) > 4:
        expcthsm = (hex2dec(range(hdrsumhex(length(hdrsumhex)-3,length(hdrsumhex))+1)))
    else:
        expcthsm = hex2dec(hdrsumhex)+1
    
    if(float(HdrCkSm)-(expcthsm)) == 0:
        
        MsgCount += 1
        tfid.write('%s\n',['Message Count = ' + str(MsgCount) + ' >> Message Type: ' + str(MssgID) + ' >>> Data Word Count: ' + str(WrdCnt)])
        #%TEST FOR A TRUNCATED MESSAGE WHILE IN HEADER OR DATA THEN. IF IT IS, THEN IT WILL
        #%NOT PARSE THE MESSAGE DATA.
        #%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
        if (atByte + 12) >= sp.nbytes:
            TruncEnd = 1
        else:
            TruncEnd = ((atByte+10+int64(WrdCnt)*2+2) > sp.nbytes)
        
        #% RUN DATA CHECKSUM FOR EVERY IDENTIFIED MESSAGE IF THERE IS
        #% CONTENT IN MESSAGE
        #%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
        if (WrdCnt != 0) and !TruncEnd:
            DataByte = atByte + 10
            MsgData = numpy.zeros(1,WrdCnt)
            for iw in range(1, int64(WrdCnt)):
                gfid.seek((DataByte + 2) * (iw-1),0)
                MsgData[iw] = fread(gfid,1, *uint16)
            gfid.seek((DataByte + 2) * iw, 0)
            DataCkSm = fread(gfid,1, *uint16)
            DataCkSm = float(DataCkSm)
            
            datsumbin = dec2bin((float(sum(MsgData))))  

            if len(datsumbin) > 16:
                if bin2dec(range(datsumbin(len(datsumbin)-15,len(datsumbin)))) != 0:
                    expctdsm = (bin2dec(range(datsumbin(len(datsumbin)-15, len(datsumbin)))))
                else:
                    expctdsm = (bin2dec(range(datsumbin(len(datsumbin)-15, len(datsumbin)))) + 2^16)
            else:
                expctdsm = bin2dec(datsumbin)
                if datsumbin !=0:
                    expctdsm += 1

            if (float(DataCkSm)+ expctdsm) == 2^16:
                tfid.write('%s\n','     Data CheckSum OK')
            else:
                DatCkSmFail += 1
                tfid.write('%s\n',' *** ERROR *** Data CheckSum FAILS!')
            #clear DataByte MsgData iw datsumbin expctdsm
            #I HAVE NO IDEA HOW TO DO THIS
        
        #% IDENTIFY MESSAGE AND APPROPRIATELY ADD TO ACCOUNTING
        #%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
        if MssgID == 3:
            Msg3Count += 1
            tfid.write('%s\n',['Mssg 3 Count: ' + str(Msg3Count)])
            gfid.seek(atByte+10, 0)
            #GPSd1(range(1,4)) = fread(gfid,4, *uint16, ieee-le)   #% GPStime
            #need to figure out how to use fread in python because this takes 4
            #arguments when I only have a predefined function for 3 arguments
            gtW1bP = dec2bin(GPSd1(2))
            gtW2bP = dec2bin(GPSd1(1))
            gtW3bP = dec2bin(GPSd1(4))
            gtW4bP = dec2bin(GPSd1(3))
            #still don't know if dec2bin exists in python
            
            if len(gtW1bP) < 16:
                n0 = 16 - len(gtW1bP)
                gtW1b = [numpy.matlib.repmat('0',1,n0), gtW1bP]
            else:
                gtW1b = gtW1bP
            
            if len(gtW2bP) < 16:
                n0 = 16 - len(gtW2bP)
                gtW2b = [numpy.matlib.repmat('0',1,n0), gtW2bP]
            else:
                gtW2b = gtW2bP
            
            if len(gtW3bP) < 16:
                n0 = 16 - len(gtW3bP)
                gtW3b = [numpy.matlib.repmat('0',1,n0), gtW3bP]
            else:
                gtW3b = gtW3bP
            
            if len(gtW4bP) < 16:
                n0 = 16 - len(gtW4bP)
                gtW4b = [numpy.matlib.repmat('0',1,n0), gtW4bP]
            else:
                gtW4b = gtW4bP

            #clear gtW1bP gtW2bP gtW3bP gtW4bP;
            gtb = [gtW1b, gtW2b, gtW3b, gtW4b]

            #%if(length(gtb) ~= 64), 'WARNING: GPS Time in binary not 64-bit!', end;

            LHSb = gtb[range(2,21)]

            LHSf = 0
            for i in range(1, len(LHSb)):
                LHSf += bin2dec(LHSb(i))*2^(-i)
            LHSf *= (2^20)

            RHSb = gtb[range(22,64)]
            RHSf = bin2dec(RHSb) * (2^(20-63))
            
            #is GPSTimeSec an array?

            GPSTimSec[Msg3Count] = LHSf + RHSf
            #clear gtW1b gtW2b gtW3b gtW4b GPSd1 LHSf RHSf;

            MsgPCntrs[1,pidx] = Msg3500PartCount
            MsgPCntrs[2,pidx] = Msg3501PartCount
            MsgPCntrs[3,pidx] = Msg3502PartCount
            MsgPCntrs[4,pidx] = Msg3512PartCount
            MsgPCntrs[5,pidx] = Msg3623PartCount
            
            Msg3500PartCount = int64(0)
            Msg3501PartCount = int64(0)
            Msg3502PartCount = int64(0)
            Msg3512PartCount = int64(0)
            Msg3623PartCount = int64(0)
            
            pidx += 1

        #%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
        if MssgID == 3500:
            Msg3500Count += 1 #%Count the number of Message 3501.
            Msg3500PartCount += 1 #%Count the number of Message 3501.
            gfid.seek(atByte+18,0)
            CurrMode3500 = fread(gfid,[1,1], int16,ieee-le) #% GPS WGS-84 latitude (semicircles)
            Msg3500CurrMode[Msg3500Count] = CurrMode3500
            tfid.write('%s\n',['Mssg 3500 Count: ' + str(Msg3Count)])
            tfid.write('%s\n',['   Current Mode: ' + cell2mat(GPSModes(CurrMode3500))]);
            #I HAVE NO IDEA HOW TO DO THIS
            #WHAT IS CELL2MAT IN PYTHON?
            #clear CurrMode3500;
            
        
        #%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
        if MssgID == 3501:
            Msg3501Count += 1 #%Count the number of Message 3501.
            Msg3501PartCount += 1 #%Count the number of Message 3501.
            
            if(PAR3501 and parsed3501==0)and not TruncEnd:
                #%D = zeros(14,1);
                
                #fix all the fread
                
                gfid.seek(atByte+10, 0)
                Gtt = fread(gfid,[4,1], 4*uint16, ieee-le) #% GPS time tag (s)
                gfid.seek(atByte+18, 0)
                Glat  = fread(gfid,[1,1], 'int32','ieee-le') #% GPS WGS-84 latitude (semicircles)
                gfid.seek(atByte+22, 0)
                Glon  = fread(gfid,[1,1], 'int32','ieee-le') #% GPS WGS-84 longitude (semicircles)
                gfid.seek(atByte+26,0)
                Galt  = fread(gfid,[1,1], 'int32','ieee-le') #% GPS WGS-84 altitude (m)
                gfid.seek(atByte+30,0)
                Gvn   = fread(gfid,[1,1], 'int32','ieee-le') #% GPS velocity N (m/s)
                gfid.seek(atByte+34, 0)
                Gve   = fread(gfid,[1,1], 'int32','ieee-le') #% GPS velocity E (m/s)
                gfid.seek(atByte+38, 0)  
                Gvu   = fread(gfid,[1,1], 'int32','ieee-le') #% GPS velocity up (m/s)
                gfid.seek(atByte+42,0)  
                Gptch = fread(gfid,[1,1], 'int32','ieee-le') #% GPS pitch (+up, semicircles)
                gfid.seek(atByte+46, 0)  
                Groll = fread(gfid,[1,1], 'int32','ieee-le') #% GPS roll (+right, semicircles)
                gfid.seek(atByte+50, 0)  
                Ghtru = fread(gfid,[1,1], 'int32','ieee-le') #% GPS true heading (+clockwise from N, semicircles))
                gfid.seek(atByte+54,0)  
                Gdchk = fread(gfid,[1,1], 'uint16'  ,'ieee-le') #% Data CheckSum
                
                GttW1bP = dec2bin(Gtt(2))
                GttW2bP = dec2bin(Gtt(1))
                GttW3bP = dec2bin(Gtt(4))
                GttW4bP = dec2bin(Gtt(3))
                
                if len(GttW1bP) < 16:
                    n0 = 16 - len(GttW1bP)
                    GttW1b = [numpy.matlib.repmat('0',1,n0), GttW1bP]
                else:
                    GttW1b = GttW1bP
                
                if len(GttW2bP) < 16:
                    n0 = 16 - len(GttW2bP)
                    GttW2b = [numpy.matlib.repmat('0',1,n0), GttW2bP]
                else:
                    GttW2b = GttW2bP
                
                if len(GttW3bP) < 16:
                    n0 = 16 - len(GttW3bP)
                    GttW3b = [numpy.matlib.repmat('0',1,n0), GttW3bP]
                else:
                    GttW3b = GttW3bP
                
                if len(GttW4bP) < 16:
                    n0 = 16 - len(GttW4bP)
                    GttW4b = [numpy.matlib.repmat('0',1,n0), GttW4bP]
                else:
                    GttW4b = GttW4bP;
                
                #clear GttW1bP GttW2bP GttW3bP GttW4bP;
                Gttb = [GttW1b, GttW2b, GttW3b, GttW4b] #around 370 in real code
                
                GttLHSb = Gttb[range(2,21)]
                GttLHSf = 0
                
                for i in range(1, len(GttLHSb)):
                     GttLHSf += bin2dec(GttLHSb(i))*2^(-i)
                
                GttLHSf *= (2^20)
                
                GttRHSb = Gttb[range(22,64)]
                GttRHSf = bin2dec(GttRHSb)*(2^(20-63))
                
                Gttsec = GttLHSf + GttRHSf

                Glat  = Glat * (180.0/(2^31))
                Glon  = Glon * (180.0/(2^31))
                Galt  = Galt * (2^(15-31))
                Gvn   = Gvn * (2^(10-31))
                Gve   = Gve * (2^(10-31))
                Gvu   = Gvu * (2^(10-31))
                Gptch = Gptch * (180.0/(2^31))
                Groll = Groll * (180.0/(2^31))
                Ghtru = Ghtru * (180.0/(2^31))

                D = [double(Msg3501Count), double(atByte), double(FlgWrd), Gttsec, Glat,\
                    Glon, Galt, Gvn, Gve, Gvu, Gptch, Groll, Ghtru, float(Gdchk)]

                fid3501.write(['%ld, %ld, %d, ' + numpy.matlib.repmat('%f, ',1,10) +  '%ld \n'],D)

                #clear D Glat Glon Galt Gvn Gve Gvu Gptch Groll Ghtru;
                #clear GttW1b GttW2b GttW3b GttW4b Gttb GttLHSb GttLHSf;
                #clear GttRHSb GttRHSf Gttsec Gdchk;
                #%keyboard
        
        
        #%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
        
        if MssgID == 3502:
            Msg3502Count += 1 #%Count the number of Message 3501.
            Msg3502PartCount += 1 #%Count the number of Message 3501.
        
        #%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
        if (MssgID == 3512) and not TruncEnd:
            Msg3512Count += 1 #%Count the number of Message 3512.
            Msg3512PartCount += 1 #%Count the number of Message 3512.
            gfid.seek(atByte+10, 0)
            GPSd2[range(1,4)] = fread(gfid,4, *uint16, ieee-le)   #% GPStime
            gtW1bP = dec2bin(GPSd2(2))
            gtW2bP = dec2bin(GPSd2(1))
            gtW3bP = dec2bin(GPSd2(4))
            gtW4bP = dec2bin(GPSd2(3))            

            if len(gtW1bP) < 16:
                n0 = 16 - len(gtW1bP)
                gtW1b = [numpy.matlib.repmat('0',1,n0), gtW1bP]
            else:
                gtW1b = gtW1bP
            
            if len(gtW2bP) < 16:
                n0 = 16 - len(gtW2bP)
                gtW2b = [numpy.matlib.repmat('0',1,n0), gtW2bP]
            else:
                gtW2b = gtW2bP
            
            if len(gtW3bP) < 16:
                n0 = 16 - len(gtW3bP)
                gtW3b = [numpy.matlib.repmat('0',1,n0), gtW3bP]
            else:
                gtW3b = gtW3bP
            
            if len(gtW4bP) < 16:
                n0 = 16 - len(gtW4bP)
                gtW4b = [numpy.matlib.repmat('0',1,n0), gtW4bP]
            else:
                gtW4b = gtW4bP
            
            #clear gtW1bP gtW2bP gtW3bP gtW4bP;
            gtb = [gtW1b, gtW2b, gtW3b, gtW4b]

            #%if(length(gtb) ~= 64), 'WARNING: GPS Time in binary not 64-bit!', end;

            LHSb = gtb[range(2,21)]

            LHSf = 0
            for i in range(1, len(LHSb)):
                LHSf += bin2dec(LHSb(i))*2^(-i)
            LHSf += (2^20)

            RHSb = gtb(range(22,64))
            RHSf = bin2dec(RHSb)*(2^(20-63))

            GPSTimSec3512[Msg3512Count] = LHSf + RHSf
            #clear gtW1b gtW2b gtW3b gtW4b GPSd2 LHSf RHSf;
            
        
        #%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
        if MssgID == 3623:
            Msg3623Count += 1 #%Count the number of Message 3501.
            Msg3623PartCount += 1 #%Count the number of Message 3501.

        #%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
        if WrdCnt > 0:
            #%if(WrdCnt > 128), '   *** WARNING *** Message Contains more than 128 words', end;
            atByte = atByte + 10 + int64(WrdCnt*2) +2        
        else:
            atByte += 10
        
        if not TruncEnd:
            tfid.write('%s \n\n',['Byte at end = ' + str(atByte)])
        else:
            tfid.write('%s \n\n',['*** CAUTION *** Last Message - Truncated'])
            atByte = sp.nbytes
        
        gfid.seek(atByte, 0) 
        WordTest = fread(gfid,1, '*uint16')
    
        while WordTest != hex2dec('81FF'):
            GarBytes += 1
            tfid.write('%s \n\n',['*** ERROR *** Excess Msg Input at Byte ' + str(atByte)])        
            atByte += 1
            gfid.seek(atByte, 0)  
            WordTest  = fread(gfid,1, '*uint16')        
        
    else: 
        GarBytes += 1
        gfid.seek(atByte+1,0) 
        WordTest = fread(gfid,1, *uint16)
    
        while WordTest != hex2dec('81FF'):
            tfid.write('%s \n\n',['*** ERROR *** False Msg Header at Byte ' + str(atByte)])
            atByte += 1
            gfid.seek(atByte, 0)  
            WordTest  = fread(gfid,1, *uint16)
            GarBytes += 1
            
            
tfid.write('%s \n','=======================END======OF======FILE==========================')

MsgPCntrs[1,pidx] = Msg3500PartCount
MsgPCntrs[2,pidx] = Msg3501PartCount
MsgPCntrs[3,pidx] = Msg3502PartCount
MsgPCntrs[4,pidx] = Msg3512PartCount
MsgPCntrs[5,pidx] = Msg3623PartCount


#%keyboard
'''
%-------------------------------------------------------------------------------------
% DETERMINE GPS MODES PRESENT AND HOW MUCH OF FILE IS IN AIR-NAVIGATION
%-------------------------------------------------------------------------------------
'''
ModeTypes  = numpy.unique(Msg3500CurrMode)
nModes = len(ModeTypes)
ModeString = cell2mat(GPSModes(ModeTypes)) #I HAVE NO IDEA HOW TO DO THIS
NavPerCent = 100.0 * len(numpy.nonzero(Msg3500CurrMode==7)) / Msg3500Count        

'''
%-------------------------------------------------------------------------------------
% CREATE THE PASS/FAIL FLAGS FOR THE THREE GPS PARAMETERS
%-------------------------------------------------------------------------------------
'''
if GarBytes == 0:
    GPSMssgPF = 'PASS'
else:
    GPSMssgPF = 'FAIL'

if DatCkSmFail == 0:
    GPSCkSmPF = 'PASS'
else:
    GPSCkSmPF = 'FAIL'
if NavPerCent == 100:
    GPSArNvPF = 'PASS'
else:
    GPSArNvPF = 'FAIL'

'''
%-------------------------------------------------------------------------------------
% POPULATE THE MSGCOUNT STRUCTURE WITH THE FINAL COUNTS
%-------------------------------------------------------------------------------------
MSGCOUNT = struct('M0003',[Msg3Count],'M3500',[Msg3500Count],'M3501',[Msg3501Count],...
                  'M3502',[Msg3502Count],'M3512',[Msg3512Count],'M3623',[Msg3623Count]);

%-------------------------------------------------------------------------------------
% POPULATE THE REPCOUNT STRUCTURE WITH THE FINAL COUNTS
%-------------------------------------------------------------------------------------
REPCOUNT = struct('TotalBytes',[atByte],'UnassignedBytes',[GarBytes],'FailedChkSum',[DatCkSmFail],...
                  'NoGPSModes',[nModes],'GPSModes',[ModeString],'PctNavMode',[NavPerCent]);
'''
'''
%-------------------------------------------------------------------------------------
% WRITE SUMMARY REPORT
%-------------------------------------------------------------------------------------
'''
ofid.write(['>>> GPS Msg. Integrity : ' + GPSMssgPF + '\n'])
ofid.write(['>>> GPS Data Integrity : ' + GPSCkSmPF + '\n'])
ofid.write(['>>> GPS AirNav Mode    : ' + GPSArNvPF + '\n'])
ofid.write('======================================================================\n')
ofid.write('SUMMARY\n')
ofid.write('======================================================================\n')
ofid.write(['Total Number of Messages 0003: ' + str(Msg3Count) + '\n'])
ofid.write(['Total Number of Messages 3500: ' + str(Msg3500Count) + '\n'])
ofid.write(['Total Number of Messages 3501: ' + str(Msg3501Count) + '\n'])
ofid.write(['Total Number of Messages 3502: ' + str(Msg3502Count) + '\n'])
ofid.write(['Total Number of Messages 3512: ' + str(Msg3512Count) + '\n'])
ofid.write(['Total Number of Messages 3623: ' + str(Msg3623Count) + '\n\n'])
ofid.write(['Bytes in file               : ' + str(atByte) + '\n'])
ofid.write(['Bytes unassigned to messages: ' + str(GarBytes) + '\n\n'])
ofid.write(['Failed Data Checksums       : ' + str(DatCkSmFail) + '\n\n'])
ofid.write(['Number of GPS Modes : ' + str(nModes) + '\n'])
ofid.write(['GPS Modes           : ' + ModeString + '\n'])
ofid.write(['Portion in Nav Mode : ' + str(round(NavPerCent)) + '%% \n'])
ofid.write('======================================================================\n')
tfid.close()
ofid.close()
'''
%-------------------------------------------------------------------------------------
% MERGE PARSED-MESSAGE FILE WITH THE SUMMARY-FILE
%-------------------------------------------------------------------------------------
'''
if ispc:
    system(['type ' + os.path.join(datadir,summfile) + ' ' + os.path.join(datadir,parsfile) + ' > ' + os.path.join(datadir,allfile)])
if isunix:
    system(['cat ' + os.path.join(datadir,summfile) + ' ' + os.path.join(datadir,parsfile) + ' > ' + os.path.join(datadir,allfile)])


os.remove(os.path.join(datadir,summfile))
os.remove(os.path.join(datadir,parsfile))