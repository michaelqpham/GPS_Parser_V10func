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

function [MSGCOUNT, REPCOUNT, GPSTimSec] = GPS_Parser_V10func(filename, datadir, GPSFigs, PAR3501)

%=====================================================================================
% U S E R    P A R A M E T E R S
%=====================================================================================
Vnum = 'V10';
PssblMsg = [3 3500 3501 3502 3512 3623];
GPSModes = {'(1) Test ','(2) Initialization ','(3) *not used* ',...
            '(4) Fine Alignment ','(5) Air Alignment ','(6) Transfer Alignment ',...
            '(7) Air Navigation ','(8) *not used* ','(9) GPS Only '};
MSGCOUNT = struct('M0003',[],'M3500',[],'M3501',[],'M3502',[],'M3512',[],'M3623',[]);
REPCOUNT = struct('TotalBytes',[],'UnassignedBytes',[],'FailedChkSum',[],...
                  'NoGPSModes',[],'GPSModes',[],'PctNavMode',[]);

%=====================================================================================

%-------------------------------------------------------------------------------------
% INITIAL USER INTERACTIONS - CHOOSE FILE AND OPT FOR FIGURE
%-------------------------------------------------------------------------------------
hgps = waitbar(0,['Processing GPS data...']);


%-------------------------------------------------------------------------------------
% OPEN THE GPS FILE
%-------------------------------------------------------------------------------------
gpsname = filename;
gfid=fopen(fullfile(datadir,gpsname), 'rb');
%cd (datadir);


%-------------------------------------------------------------------------------------
% DETERMINE IF MSG3501 ALREADY EXISTS AND IF NOT, CREATE IT.
%-------------------------------------------------------------------------------------
msg3501Filename = strcat(strrep(gpsname,'.bin',''),'-Msg3501.csv');
parsed3501 = exist(fullfile(datadir,msg3501Filename),'file');    

if ((PAR3501==true)&&(parsed3501==0))
    fid3501=fopen(fullfile(datadir,msg3501Filename),'w+');
    fprintf(fid3501,'%s \n', '%=============================================================================');
    fprintf(fid3501,'%s \n',['%NGDCS Check ' Vnum ': Message 3501 Navigation Solution of ' gpsname]);
    fprintf(fid3501,'%s \n', '%* If present "Date" and "Abs Time" are derived and not original Msg3501 data');
    fprintf(fid3501,'%s \n', '%=============================================================================');
    fprintf(fid3501,'%s \n', '%Msg3501#, StartByte, Flag, Time, Lat, Lon, Alt, VelN, VelE, VelUp, Pitch, Roll, HeadTrue, DataCheckSum');
end

%-------------------------------------------------------------------------------------
% READ C-MIGITS III DATA AND FIND MESSAGE ID3 FOR EACH PULSE. BUT FIRST,
% SEEK THE FIRST hx81FF CODE FOR THE FIRST COMPLETE MESSAGE.
%-------------------------------------------------------------------------------------
allfile  = strcat(filename,'-parsed.txt');
summfile = strcat(filename,'-summtemp.txt');
parsfile = strcat(filename,'-parsetemp.txt');

ofid=fopen(fullfile(datadir,summfile),'w+');
tfid=fopen(fullfile(datadir,parsfile),'w+');

sp = dir(fullfile(datadir,gpsname));
nwrd = sp.bytes/2;
nams = int64(nwrd/5);

GPSTimSec = zeros(1,nams);
GPSTimSec3512 = zeros(1,nams);
MsgPCntrs = zeros(5,nams);

Estampa = ['GPS_Parser_' Vnum '  on ' datestr(now)];
fprintf(ofid,'%s \n','======================================================================');
fprintf(ofid,'%s \n\n', ['   ' Estampa]);
fprintf(ofid,'%s \n',['   GPS File: ' gpsname]);
fprintf(ofid,'%s \n','======================================================================');

atByte = int64(0);
MsgCount = int64(0);
Msg3Count = int64(0);
Msg3500Count = int64(0);
Msg3501Count = int64(0);
Msg3502Count = int64(0);
Msg3512Count = int64(0);
Msg3623Count = int64(0);
Msg3500PartCount = int64(0);
Msg3501PartCount = int64(0);
Msg3502PartCount = int64(0);
Msg3512PartCount = int64(0);
Msg3623PartCount = int64(0);
Msg3500CurrMode  = int64(0);
pidx = 1;
DatCkSmFail = 0;

fseek(gfid,0,'bof');
WordTest = (fread(gfid,1, '*uint16'));
while (WordTest ~= hex2dec('81FF'))
    atByte = atByte + 1;
    fseek(gfid,atByte,'bof');  WordTest  = fread(gfid,1, '*uint16');
    %WordTest = dec2hex(WordTest)
end

GarBytes = int64(0);
count = int64(0);
while (~feof(gfid))
    atByte;
    if (atByte>count+100000)
        count = atByte;
        %pause(1);
        waitbar(double(count)/sp.bytes,hgps,['Processing GPS data...' num2str(100*count/sp.bytes) '%']);
    end
    fprintf(tfid,'%s \n','------------------------------');
    fprintf(tfid,'%s \n',['Byte at start = ' num2str(atByte)]);
    
    fseek(gfid,atByte+0,'bof');  SynWrd  = fread(gfid,1, '*uint16');
    fseek(gfid,atByte+2,'bof');  MssgID  = fread(gfid,1, '*uint16');             % Message 3 ID
    fseek(gfid,atByte+4,'bof');  WrdCnt  = fread(gfid,1, '*uint16');             % Word count
    fseek(gfid,atByte+6,'bof');  FlgWrd  = fread(gfid,1, '*uint16');             % Flag word
    fseek(gfid,atByte+8,'bof');  HdrCkSm = fread(gfid,1, '*uint16');             % Header Check Sum
    
    hdrsumhex = dec2hex(bitcmp(uint32(SynWrd)+uint32(MssgID)+uint32(WrdCnt)+uint32(FlgWrd)));
    if(length(hdrsumhex)>4)
        expcthsm = (hex2dec(hdrsumhex(length(hdrsumhex)-3:length(hdrsumhex)))+1);
    else
        expcthsm = hex2dec(hdrsumhex)+1;
    end
    
    if((double(HdrCkSm)-(expcthsm))==0)
        
        MsgCount = MsgCount + 1;
        fprintf(tfid,'%s\n',['Message Count = ' num2str(MsgCount) ' >> Message Type: ' num2str(MssgID) ' >>> Data Word Count: ' num2str(WrdCnt)]);

        %TEST FOR A TRUNCATED MESSAGE WHILE IN HEADER OR DATA THEN. IF IT IS, THEN IT WILL
        %NOT PARSE THE MESSAGE DATA.
        %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
        if(atByte+12 >= sp.bytes)
            TruncEnd = 1;
        else
            TruncEnd = ((atByte+10+int64(WrdCnt)*2+2)>sp.bytes);
        end
        
        % RUN DATA CHECKSUM FOR EVERY IDENTIFIED MESSAGE IF THERE IS
        % CONTENT IN MESSAGE
        %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
        if((WrdCnt ~= 0)&&(~TruncEnd))
            DataByte = atByte+10;
            MsgData = (zeros(1,WrdCnt));
            for iw = 1:int64(WrdCnt)
                fseek(gfid,DataByte+2*(iw-1),'bof'); MsgData(iw) = fread(gfid,1, '*uint16');
            end
            fseek(gfid,DataByte+2*iw,'bof'); DataCkSm = fread(gfid,1, '*uint16');
            DataCkSm = double(DataCkSm);
            
            datsumbin = dec2bin((double(sum(MsgData))));  

            if(length(datsumbin)>16)
                if(bin2dec((datsumbin(length(datsumbin)-15:length(datsumbin))))~=0)
                    expctdsm = (bin2dec(datsumbin(length(datsumbin)-15:length(datsumbin))));
                else
                    expctdsm = (bin2dec(datsumbin(length(datsumbin)-15:length(datsumbin)))) + 2^16;
                end                
            else
                expctdsm = bin2dec(datsumbin);
                if (datsumbin~=0)
                    expctdsm = expctdsm + 1;
                end
            end

            if((double(DataCkSm)+(expctdsm)) == 2^16)
                fprintf(tfid,'%s\n','     Data CheckSum OK');
            else
                DatCkSmFail = DatCkSmFail + 1;
                fprintf(tfid,'%s\n',' *** ERROR *** Data CheckSum FAILS!');
            end
            clear DataByte MsgData iw datsumbin expctdsm;
        end
        
        % IDENTIFY MESSAGE AND APPROPRIATELY ADD TO ACCOUNTING
        %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
        if (MssgID == 3)
            Msg3Count = Msg3Count + 1;
            fprintf(tfid,'%s\n',['Mssg 3 Count: ' num2str(Msg3Count)]);
            fseek(gfid,atByte+10,'bof');   GPSd1(1:4) = fread(gfid,4, '*uint16','ieee-le');   % GPStime
            gtW1bP = dec2bin(GPSd1(2));
            gtW2bP = dec2bin(GPSd1(1));
            gtW3bP = dec2bin(GPSd1(4));
            gtW4bP = dec2bin(GPSd1(3));
            
            if(length(gtW1bP) < 16);
                n0 = 16-length(gtW1bP);
                gtW1b = [repmat('0',1,n0) gtW1bP];
            else
                gtW1b = gtW1bP;
            end
            
            if(length(gtW2bP) < 16);
                n0 = 16-length(gtW2bP);
                gtW2b = [repmat('0',1,n0) gtW2bP];
            else
                gtW2b = gtW2bP;
            end
            
            if(length(gtW3bP) < 16);
                n0 = 16-length(gtW3bP);
                gtW3b = [repmat('0',1,n0) gtW3bP];
            else
                gtW3b = gtW3bP;
            end
            
            if(length(gtW4bP) < 16);
                n0 = 16-length(gtW4bP);
                gtW4b = [repmat('0',1,n0) gtW4bP];
            else
                gtW4b = gtW4bP;
            end

            clear gtW1bP gtW2bP gtW3bP gtW4bP;
            gtb = [gtW1b gtW2b gtW3b gtW4b];

            %if(length(gtb) ~= 64), 'WARNING: GPS Time in binary not 64-bit!', end;

            LHSb = gtb(2:21);

            LHSf = 0;
            for i = 1:length(LHSb)
                LHSf = LHSf + bin2dec(LHSb(i))*2^(-i);
            end
            LHSf = LHSf*(2^20);

            RHSb = gtb(22:64);
            RHSf = bin2dec(RHSb)*(2^(20-63));

            GPSTimSec(Msg3Count) = LHSf+RHSf;
            clear gtW1b gtW2b gtW3b gtW4b GPSd1 LHSf RHSf;

            MsgPCntrs(1,pidx) = Msg3500PartCount;
            MsgPCntrs(2,pidx) = Msg3501PartCount;
            MsgPCntrs(3,pidx) = Msg3502PartCount;
            MsgPCntrs(4,pidx) = Msg3512PartCount;
            MsgPCntrs(5,pidx) = Msg3623PartCount;
            
            Msg3500PartCount = int64(0);
            Msg3501PartCount = int64(0);
            Msg3502PartCount = int64(0);
            Msg3512PartCount = int64(0);
            Msg3623PartCount = int64(0);
            
            pidx = pidx + 1;
        end

        %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
        if (MssgID == 3500)
            Msg3500Count = Msg3500Count + 1; %Count the number of Message 3501.
            Msg3500PartCount = Msg3500PartCount + 1; %Count the number of Message 3501.
            fseek(gfid,atByte+18,'bof'); CurrMode3500  = fread(gfid,[1,1], 'int16','ieee-le'); % GPS WGS-84 latitude (semicircles)
            Msg3500CurrMode(Msg3500Count) = CurrMode3500;
            fprintf(tfid,'%s\n',['Mssg 3500 Count: ' num2str(Msg3Count)]);
            fprintf(tfid,'%s\n',['   Current Mode: ' cell2mat(GPSModes(CurrMode3500))]);
            clear CurrMode3500;
            
        end
        
        %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
        if (MssgID == 3501)
            Msg3501Count = Msg3501Count + 1; %Count the number of Message 3501.
            Msg3501PartCount = Msg3501PartCount + 1; %Count the number of Message 3501.
            
            if((PAR3501==true)&&(parsed3501==0)&&(~TruncEnd))
                %D = zeros(14,1);
                fseek(gfid,atByte+10,'bof');  Gtt   = fread(gfid,[4,1], '4*uint16','ieee-le'); % GPS time tag (s)
                fseek(gfid,atByte+18,'bof');  Glat  = fread(gfid,[1,1], 'int32','ieee-le'); % GPS WGS-84 latitude (semicircles)
                fseek(gfid,atByte+22,'bof');  Glon  = fread(gfid,[1,1], 'int32','ieee-le'); % GPS WGS-84 longitude (semicircles)
                fseek(gfid,atByte+26,'bof');  Galt  = fread(gfid,[1,1], 'int32','ieee-le'); % GPS WGS-84 altitude (m)
                fseek(gfid,atByte+30,'bof');  Gvn   = fread(gfid,[1,1], 'int32','ieee-le'); % GPS velocity N (m/s)
                fseek(gfid,atByte+34,'bof');  Gve   = fread(gfid,[1,1], 'int32','ieee-le'); % GPS velocity E (m/s)
                fseek(gfid,atByte+38,'bof');  Gvu   = fread(gfid,[1,1], 'int32','ieee-le'); % GPS velocity up (m/s)
                fseek(gfid,atByte+42,'bof');  Gptch = fread(gfid,[1,1], 'int32','ieee-le'); % GPS pitch (+up, semicircles)
                fseek(gfid,atByte+46,'bof');  Groll = fread(gfid,[1,1], 'int32','ieee-le'); % GPS roll (+right, semicircles)
                fseek(gfid,atByte+50,'bof');  Ghtru = fread(gfid,[1,1], 'int32','ieee-le'); % GPS true heading (+clockwise from N, semicircles))
                fseek(gfid,atByte+54,'bof');  Gdchk = fread(gfid,[1,1], 'uint16'  ,'ieee-le'); % Data CheckSum
                
                GttW1bP = dec2bin(Gtt(2));
                GttW2bP = dec2bin(Gtt(1));
                GttW3bP = dec2bin(Gtt(4));
                GttW4bP = dec2bin(Gtt(3));
                
                if(length(GttW1bP) < 16)
                    n0 = 16-length(GttW1bP);
                    GttW1b = [repmat('0',1,n0) GttW1bP];
                else
                    GttW1b = GttW1bP;
                end
                
                if(length(GttW2bP) < 16)
                    n0 = 16-length(GttW2bP);
                    GttW2b = [repmat('0',1,n0) GttW2bP];
                else
                    GttW2b = GttW2bP;
                end
                
                if(length(GttW3bP) < 16)
                    n0 = 16-length(GttW3bP);
                    GttW3b = [repmat('0',1,n0) GttW3bP];
                else
                    GttW3b = GttW3bP;
                end
                
                if(length(GttW4bP) < 16)
                    n0 = 16-length(GttW4bP);
                    GttW4b = [repmat('0',1,n0) GttW4bP];
                else
                    GttW4b = GttW4bP;
                end
                
                clear GttW1bP GttW2bP GttW3bP GttW4bP;
                Gttb = [GttW1b GttW2b GttW3b GttW4b];
                
                GttLHSb = Gttb(2:21);
                GttLHSf = 0;
                
                for i = 1:length(GttLHSb)
                     GttLHSf = GttLHSf + bin2dec(GttLHSb(i))*2^(-i);
                end
                
                GttLHSf = GttLHSf*(2^20);
                
                GttRHSb = Gttb(22:64);
                GttRHSf = bin2dec(GttRHSb)*(2^(20-63));
                
                Gttsec = GttLHSf+GttRHSf;

                Glat  = Glat*(180/(2^31));
                Glon  = Glon*(180/(2^31));
                Galt  = Galt*(2^(15-31));
                Gvn   = Gvn*(2^(10-31));
                Gve   = Gve*(2^(10-31));
                Gvu   = Gvu*(2^(10-31));
                Gptch = Gptch*(180/(2^31));
                Groll = Groll*(180/(2^31));
                Ghtru = Ghtru*(180/(2^31));

                D = [double(Msg3501Count); double(atByte); double(FlgWrd); Gttsec; Glat;...
                    Glon; Galt; Gvn; Gve; Gvu; Gptch; Groll; Ghtru; double(Gdchk)];

                fprintf(fid3501,['%ld, %ld, %d, ' repmat('%f, ',1,10) '%ld \n'],D);

                clear D Glat Glon Galt Gvn Gve Gvu Gptch Groll Ghtru;
                clear GttW1b GttW2b GttW3b GttW4b Gttb GttLHSb GttLHSf;
                clear GttRHSb GttRHSf Gttsec Gdchk;
                %keyboard
            end
        end
        
        
        %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
        
        if (MssgID == 3502)
            Msg3502Count = Msg3502Count + 1; %Count the number of Message 3501.
            Msg3502PartCount = Msg3502PartCount + 1; %Count the number of Message 3501.
        end
        
        %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
        if ((MssgID == 3512)&&(~TruncEnd))
            Msg3512Count = Msg3512Count + 1; %Count the number of Message 3512.
            Msg3512PartCount = Msg3512PartCount + 1; %Count the number of Message 3512.
            fseek(gfid,atByte+10,'bof');   GPSd2(1:4) = fread(gfid,4, '*uint16','ieee-le');   % GPStime
            gtW1bP = dec2bin(GPSd2(2));
            gtW2bP = dec2bin(GPSd2(1));
            gtW3bP = dec2bin(GPSd2(4));
            gtW4bP = dec2bin(GPSd2(3));            

            if(length(gtW1bP) < 16);
                n0 = 16-length(gtW1bP);
                gtW1b = [repmat('0',1,n0) gtW1bP];
            else
                gtW1b = gtW1bP;
            end
            
            if(length(gtW2bP) < 16);
                n0 = 16-length(gtW2bP);
                gtW2b = [repmat('0',1,n0) gtW2bP];
            else
                gtW2b = gtW2bP;
            end
            
            if(length(gtW3bP) < 16);
                n0 = 16-length(gtW3bP);
                gtW3b = [repmat('0',1,n0) gtW3bP];
            else
                gtW3b = gtW3bP;
            end
            
            if(length(gtW4bP) < 16);
                n0 = 16-length(gtW4bP);
                gtW4b = [repmat('0',1,n0) gtW4bP];
            else
                gtW4b = gtW4bP;
            end
            
            clear gtW1bP gtW2bP gtW3bP gtW4bP;
            gtb = [gtW1b gtW2b gtW3b gtW4b];

            %if(length(gtb) ~= 64), 'WARNING: GPS Time in binary not 64-bit!', end;

            LHSb = gtb(2:21);

            LHSf = 0;
            for i = 1:length(LHSb)
                LHSf = LHSf + bin2dec(LHSb(i))*2^(-i);
            end
            LHSf = LHSf*(2^20);

            RHSb = gtb(22:64);
            RHSf = bin2dec(RHSb)*(2^(20-63));

            GPSTimSec3512(Msg3512Count) = LHSf+RHSf;
            clear gtW1b gtW2b gtW3b gtW4b GPSd2 LHSf RHSf;
            
        end
        
        %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
        if (MssgID == 3623)
            Msg3623Count = Msg3623Count + 1; %Count the number of Message 3501.
            Msg3623PartCount = Msg3623PartCount + 1; %Count the number of Message 3501.
        end

        %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
        if (WrdCnt > 0)
            %if(WrdCnt > 128), '   *** WARNING *** Message Contains more than 128 words', end;
            atByte = atByte + 10 + int64(WrdCnt*2) +2;        
        else
            atByte = atByte + 10;
        end
        
        if(~TruncEnd)
            fprintf(tfid,'%s \n\n',['Byte at end = ' num2str(atByte)]);
        else
            fprintf(tfid,'%s \n\n',['*** CAUTION *** Last Message - Truncated']);
            atByte = sp.bytes;
        end
        
        fseek(gfid,atByte,'bof'); WordTest = (fread(gfid,1, '*uint16'));
    
        while (WordTest ~= hex2dec('81FF'))
            GarBytes = GarBytes + 1;
            fprintf(tfid,'%s \n\n',['*** ERROR *** Excess Msg Input at Byte ' num2str(atByte)]);        
            atByte = atByte + 1;
            fseek(gfid,atByte,'bof');  WordTest  = fread(gfid,1, '*uint16');        
        end
        
    else        
        GarBytes = GarBytes + 1;
        fseek(gfid,atByte+1,'bof'); WordTest = (fread(gfid,1, '*uint16'));
    
        while (WordTest ~= hex2dec('81FF'))
            fprintf(tfid,'%s \n\n',['*** ERROR *** False Msg Header at Byte ' num2str(atByte)]);
            atByte = atByte + 1;
            fseek(gfid,atByte,'bof');  WordTest  = fread(gfid,1, '*uint16');
            GarBytes = GarBytes + 1;
        end
    end
        
end

fprintf(tfid,'%s \n','=======================END======OF======FILE==========================');

MsgPCntrs(1,pidx) = Msg3500PartCount;
MsgPCntrs(2,pidx) = Msg3501PartCount;
MsgPCntrs(3,pidx) = Msg3502PartCount;
MsgPCntrs(4,pidx) = Msg3512PartCount;
MsgPCntrs(5,pidx) = Msg3623PartCount;


%keyboard
%-------------------------------------------------------------------------------------
% DETERMINE GPS MODES PRESENT AND HOW MUCH OF FILE IS IN AIR-NAVIGATION
%-------------------------------------------------------------------------------------
ModeTypes  = unique(Msg3500CurrMode);
nModes = length(ModeTypes);
ModeString = cell2mat(GPSModes(ModeTypes));
NavPerCent = 100*length(find(Msg3500CurrMode==7))/Msg3500Count;        


%-------------------------------------------------------------------------------------
% CREATE THE PASS/FAIL FLAGS FOR THE THREE GPS PARAMETERS
%-------------------------------------------------------------------------------------

if (GarBytes==0)
    GPSMssgPF = 'PASS';
else
    GPSMssgPF = 'FAIL';
end

if (DatCkSmFail==0)
    GPSCkSmPF = 'PASS';
else
    GPSCkSmPF = 'FAIL';
end
if (NavPerCent==100)
    GPSArNvPF = 'PASS';
else
    GPSArNvPF = 'FAIL';
end


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


%-------------------------------------------------------------------------------------
% WRITE SUMMARY REPORT
%-------------------------------------------------------------------------------------
fprintf(ofid, ['>>> GPS Msg. Integrity : ' GPSMssgPF  '\n']);
fprintf(ofid, ['>>> GPS Data Integrity : ' GPSCkSmPF '\n']);
fprintf(ofid, ['>>> GPS AirNav Mode    : ' GPSArNvPF '\n']);

fprintf(ofid, '======================================================================\n');
fprintf(ofid, 'SUMMARY\n');
fprintf(ofid, '======================================================================\n');
fprintf(ofid, ['Total Number of Messages 0003: ' num2str(Msg3Count) '\n']);
fprintf(ofid, ['Total Number of Messages 3500: ' num2str(Msg3500Count) '\n']);
fprintf(ofid, ['Total Number of Messages 3501: ' num2str(Msg3501Count) '\n']);
fprintf(ofid, ['Total Number of Messages 3502: ' num2str(Msg3502Count) '\n']);
fprintf(ofid, ['Total Number of Messages 3512: ' num2str(Msg3512Count) '\n']);
fprintf(ofid, ['Total Number of Messages 3623: ' num2str(Msg3623Count) '\n\n']);
fprintf(ofid, ['Bytes in file               : ' num2str(atByte) '\n']);
fprintf(ofid, ['Bytes unassigned to messages: ' num2str(GarBytes) '\n\n']);
fprintf(ofid, ['Failed Data Checksums       : ' num2str(DatCkSmFail) '\n\n']);
fprintf(ofid, ['Number of GPS Modes : ' num2str(nModes) '\n']);
fprintf(ofid, ['GPS Modes           : ' ModeString '\n']);
fprintf(ofid, ['Portion in Nav Mode : ' num2str(round(NavPerCent)) '%% \n']);
fprintf(ofid, '======================================================================\n');

fclose(tfid);
fclose(ofid);

%-------------------------------------------------------------------------------------
% MERGE PARSED-MESSAGE FILE WITH THE SUMMARY-FILE
%-------------------------------------------------------------------------------------
if(ispc);
    system(['type ' fullfile(datadir,summfile) ' ' fullfile(datadir,parsfile) ' > ' fullfile(datadir,allfile)]);
end
if(isunix)
    system(['cat ' fullfile(datadir,summfile) ' ' fullfile(datadir,parsfile) ' > ' fullfile(datadir,allfile)]);
end


delete(fullfile(datadir,summfile));
delete(fullfile(datadir,parsfile));


%-------------------------------------------------------------------------------------
% CREATE FIGURES IF SO DESIRED
%-------------------------------------------------------------------------------------
if(GPSFigs == true)
    figname = strcat(filename,'-parsed-fig');
    
    %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
    
    figure('Position',[10 50 600 900],'PaperUnits','inches','PaperPosition',...
         [0.5,0.5,7,10]);
    subplot(2,1,1);
    plot(GPSTimSec(1:Msg3Count),'-ko','MarkerFaceColor','k','MarkerSize',3);
    xlabel('MSG ID 3 Count','FontSize',12,'FontWeight','demi');
    ylabel('GPS Time (sec)','FontSize',12,'FontWeight','demi');
    title(filename,'Interpreter','none','FontSize',14,'FontWeight','demi');
    set(gca,'PlotBoxAspectRatio',[1.5 1 1],'TickLength',[.02 .02],...
            'XMinorTick','on','YMinorTick','on','YDir','normal');
    subplot(2,1,2);
    plot(diff(GPSTimSec(1:Msg3Count)),'-ko','MarkerFaceColor','k','MarkerSize',3);
    xlabel('MSG ID 3 Count','FontSize',12,'FontWeight','demi');
    ylabel('Delta GPS Time (sec)','FontSize',12,'FontWeight','demi');
    set(gca,'PlotBoxAspectRatio',[1.5 1 1],'TickLength',[.02 .02],...
            'XMinorTick','on','YMinorTick','on','YDir','normal');
    %
    %print('-dpsc2', '-r600', '-append', figname);
    saveas(gcf,fullfile(datadir,strcat(filename,'-gpstime-fig.pdf')),'pdf');
    
    %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
    
    figure('Position',[10 50 600 900],'PaperUnits','inches','PaperPosition',...
         [0.5,0.5,7,10]);
    subplot(2,1,1);
    plot(MsgPCntrs(1,1:pidx),'-bo','MarkerFaceColor','b','MarkerSize',3);
    xlabel('MSG ID 3 Count','FontSize',12,'FontWeight','demi');
    ylabel('Number of Msg. 3500 per Msg. 3','FontSize',12,'FontWeight','demi');
    title(filename,'Interpreter','none','FontSize',14,'FontWeight','demi');
    set(gca,'PlotBoxAspectRatio',[1.5 1 1],'TickLength',[.02 .02],...
            'XMinorTick','on','YMinorTick','on','YDir','normal');
    subplot(2,1,2);
    plot(MsgPCntrs(2,1:pidx),'-ro','MarkerFaceColor','r','MarkerSize',3);
    xlabel('MSG ID 3 Count','FontSize',12,'FontWeight','demi');
    ylabel('Number of Msg. 3501 per Msg. 3','FontSize',12,'FontWeight','demi');
    title(filename,'Interpreter','none','FontSize',14,'FontWeight','demi');
    set(gca,'PlotBoxAspectRatio',[1.5 1 1],'TickLength',[.02 .02],...
            'XMinorTick','on','YMinorTick','on','YDir','normal');
    %
    saveas(gcf,fullfile(datadir,strcat(filename,'-3500.3501-fig.pdf')),'pdf');
    %saveas(gcf,[filename '-parser-fig.pdf'],'pdf');
    %print('-dpsc2', '-r600', '-append', figname);
    
    %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
    
    figure('Position',[10 50 600 900],'PaperUnits','inches','PaperPosition',...
         [0.5,0.5,7,10]);
    subplot(2,1,1);
    plot(MsgPCntrs(3,1:pidx),'-go','MarkerFaceColor','g','MarkerSize',3);
    xlabel('MSG ID 3 Count','FontSize',12,'FontWeight','demi');
    ylabel('Number of Msg. 3502 per Msg. 3','FontSize',12,'FontWeight','demi');
    title(filename,'Interpreter','none','FontSize',14,'FontWeight','demi');
    set(gca,'PlotBoxAspectRatio',[1.5 1 1],'TickLength',[.02 .02],...
            'XMinorTick','on','YMinorTick','on','YDir','normal');
    %set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[.02 .02])
    subplot(2,1,2);
    plot(MsgPCntrs(4,1:pidx),'-mo','MarkerFaceColor','m','MarkerSize',3);
    xlabel('MSG ID 3 Count','FontSize',12,'FontWeight','demi');
    ylabel('Number of Msg. 3512 per Msg. 3','FontSize',12,'FontWeight','demi');
    title(filename,'Interpreter','none','FontSize',14,'FontWeight','demi');
    set(gca,'PlotBoxAspectRatio',[1.5 1 1],'TickLength',[.02 .02],...
            'XMinorTick','on','YMinorTick','on','YDir','normal');
    saveas(gcf,fullfile(datadir,strcat(filename,'-3502.3512-fig.pdf')),'pdf');
    %saveas(gcf,[filename '-parser-fig.pdf'],'pdf');
    %print('-dpsc2', '-r600', '-append', figname);

    %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
    
    figure('Position',[10 50 600 900],'PaperUnits','inches','PaperPosition',...
         [0.5,0.5,7,10]);
    subplot(2,1,1);
    plot(diff(GPSTimSec3512(1:Msg3512Count)),'-mo','MarkerFaceColor','c','MarkerSize',3);
    xlabel('MSG ID 3512 Count','FontSize',12,'FontWeight','demi');
    ylabel('3512 GPS Time Step (sec)','FontSize',12,'FontWeight','demi');
    title(filename,'Interpreter','none','FontSize',14,'FontWeight','demi');
    set(gca,'PlotBoxAspectRatio',[1.5 1 1],'TickLength',[.02 .02],...
            'XMinorTick','on','YMinorTick','on','YDir','normal');
    subplot(2,1,2);
    plot(MsgPCntrs(1,1:pidx),'-co','MarkerFaceColor','c','MarkerSize',3);
    xlabel('MSG ID 3 Count','FontSize',12,'FontWeight','demi');
    ylabel('Number of Msg. 3623 per Msg. 3','FontSize',12,'FontWeight','demi');
    title(filename,'Interpreter','none','FontSize',14,'FontWeight','demi');
    set(gca,'PlotBoxAspectRatio',[1.5 1 1],'TickLength',[.02 .02],...
            'XMinorTick','on','YMinorTick','on','YDir','normal');
    saveas(gcf,fullfile(datadir,strcat(filename,'-3512.3623-fig.pdf')),'pdf');
    set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[.02 .02]);
    %print('-dpsc2', '-r600', '-append', figname);
    close all;%
    
    %ps2pdf('psfile',[figname '.ps'],'pdffile',...
    %        [figname '.pdf'],'gspapersize','letter','deletepsfile', 1)
    
end

fclose all;

%-------------------------------------------------------------------------------------
% FOR 3501MSG FILE, ADD DATE AND ABSOLUTE TIME COLUMNS
%-------------------------------------------------------------------------------------
if(PAR3501 == true)
    msg3501Filename = strcat(strrep(filename,'.bin',''),'-Msg3501.csv');
    generateAbsTimeFunc( datadir, msg3501Filename );
end

close(hgps);








