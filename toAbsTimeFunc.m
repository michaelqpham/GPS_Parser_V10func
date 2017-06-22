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
% Version 11
% Daniel Nunes, Didier Keymeulen, Jacky Lee, JPL/Caltech
% DCN, 2015/MAR/27
% Daniel.Nunes@jpl.nasa.gov
%
%--------------------------------------------------------------------------
%
% Matlab function returns an absolute time of the week after converting the
% relative time in decimal seconds from a the 3501 message
%
%[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]

function dayAndTime = toAbsTimeFunc( timeRelative )

secPerDay = 24 * 60 * 60;

% assign integer and string representation of day of the week
intDayOfWeek = fix(timeRelative/secPerDay);
% DO NOT NEED UNLESS YOU WANT THE RETURN STRING TO INCLUDE DAY OF THE WEEK
% switch(intDayOfWeek)
%     case 1 
%         strDay = 'Mon';
%     case 2
%         strDay = 'Tue';
%     case 3
%         strDay = 'Wed';
%     case 4
%         strDay = 'Thu';
%     case 5
%         strDay = 'Fri';
%     case 6
%         strDay = 'Sat';
%     case 7
%         strDay = 'Sun';
% end

%hours elapsed from start of the day
secElapsed = timeRelative - (intDayOfWeek * secPerDay);
hour = fix(secElapsed/(60*60));

%minute elapsed from the hour
secElapsed = secElapsed - (hour*60*60); 
minute = fix(secElapsed/(60));

%second elapsed from the minute
second = secElapsed - (minute*60);

time = [intDayOfWeek, hour, minute, second];

%dayAndTime = sprintf('%.0f%.0f%.6f',time);
dayAndTime = time;

