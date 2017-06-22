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
% Daniel Nunes, Didier Keymeulen, Jacky Lee, JPL/Caltech, Michael Pham
% DCN, 2015/FEB/26
% Daniel.Nunes@jpl.nasa.gov
%
%--------------------------------------------------------------------------
%
% Matlab script to run the GPS_Parser; bulk of the work done by calling
% GPS_Parser function
%
%[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]
clear all
close all

%==========================================================================
% U S E R    P A R A M E T E R S
%==========================================================================
Vnum = 'V10';

PssblMsg = [3 3500 3501 3502 3512 3623];
GPSModes = {'(1) Test ','(2) Initialization ','(3) *not used* ',...
            '(4) Fine Alignment ','(5) Air Alignment ','(6) Transfer Alignment ',...
            '(7) Air Navigation ','(8) *not used* ','(9) GPS Only '};
%==========================================================================

%--------------------------------------------------------------------------
% INITIAL USER INTERACTIONS - CHOOSE FILE AND OPT FOR FIGURE
%--------------------------------------------------------------------------

% select directory and list file 
[filename, datadir] = uigetfile({'*','All Files'} ,['GPS_Parser ' Vnum ' SELECT GPS or NAV FILE'] );
%cd(datadir);

% want figure?
qstring = 'Output figure with GPS time from Mssg ID 3s ?';
choice = questdlg(qstring,'FIGURE OUTPUT','Yes','No','Yes');
if (strcmp(choice,'Yes'))
    GPSFigs = true;
else
    GPSFigs = false;    
end

% extract GPS data?
qstring = 'Extract Msg3501 GPS data to CSV file?';
choice = questdlg(qstring,'FIGURE OUTPUT','Yes','No','Yes');
if (strcmp(choice,'Yes'))
    PAR3501 = true;
else
    PAR3501 = false;    
end

%function generates the figures and parse file depending on user choice as 
% well as return data needed for Check_Flight script
[MSGCOUNT, REPCOUNT, GPSTimSec] = GPS_Parser_V10func(filename, datadir, GPSFigs, PAR3501);




