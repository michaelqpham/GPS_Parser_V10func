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
% Function takes a user specified Msg3501.csv and generate a new file named
% Msg3501AbsTime.csv with the relative time covereted to absolute time
%
%[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]
function  generateAbsTimeFunc( dataDir, readFile )

% open the read and write file
readFileID = fopen(fullfile(dataDir,readFile), 'r');
writeFile = [readFile(1:end-4) 'AbsTime.csv'];
writeFileID = fopen(fullfile(dataDir,writeFile), 'w');

% writes HEADERS that are above the column data headers
% - data header by finding 'Msg3501' starting in either the first
%   or second index in each line of text
readLine = fgetl(readFileID);
while ~(strcmp(readLine(1:7),'Msg3501') || strcmp(readLine(2:8),'Msg3501') )
    fprintf(writeFileID, '%s\n', readLine);
    readLine = fgetl(readFileID);
end

% write COLUMN HEADERS
% Inserts an additional date column 
    tmp = textscan(readLine,'%s','Delimiter',','); 
    readLine =  tmp{1}'; 
    
    %Don't try to convert if there is already a date column in the 
    % 3501msg file. It is already converted
    if( strcmp(cell2mat(readLine(5)),'Date'))
        fclose('all');
        delete(fullfile(dataDir,writeFile));
        return;
    end
    
    %leave first 4 columns unchanged
    for i = 1:4
        fprintf(writeFileID, '%s, ', cell2mat(readLine(i)));
    end
    
    %add a date column as the 5th column
        fprintf(writeFileID, 'Date, ');
        
    %add a absolute time column as the 5th column
    fprintf(writeFileID, 'Abs Time, ');
    
    %leave rest of the columns unchanged
    for i = 5:13
        fprintf(writeFileID, '%s, ', cell2mat(readLine(i)));
    end
        fprintf(writeFileID, '%s\n', cell2mat(readLine(14)));

        
%Determine the date to be inserted into the enewly created date column 
    indexFound = strfind(readFile,'20'); %looks for the leading 2 digit of the year
    date = readFile(indexFound(1):indexFound+7); 
    
    %may or may not need these variable depending on how we will increment
    %the day should a flight elapse through midnight
    yyyy = str2double(date(1:4));
    mm = str2double(date(5:6));
    dd = str2double(date(7:8));
    dowi = weekday([num2str(yyyy) '-' num2str(mm) '-' num2str(dd)]);
    
%     %Determine the day of the week
%     %SUNDAY = 1, MON = 2, ... SAT = 6
%     dayNum = weekday(datenum('20140613','yyyymmdd')); %**** I AM GETTING WRONG DAY OF WEEK
% 
%     %Convert dayNum to MON = 1, TUE = 2,... SUN = 7
%     if dayNum == 1
%         dayNum = 7;
%     else
%         dayNum = dayNum - 1;
%     end
     
% Write REMAINING ROWS of data
% -inserts a date into the date column  
% -convert relative time to absolute time
while ~feof(readFileID)
    tmp = textscan(fgetl(readFileID),'%s','Delimiter',',');
    readLine =  tmp{1}';
    
    %leave first 4 columns unchanged
    for i = 1:4
        fprintf(writeFileID, '%s, ', cell2mat(readLine(i)));
    end
    
    %convert the relative time to absolute time
    absTime = toAbsTimeFunc(str2double(cell2mat(readLine(4))));
    dow = absTime(1);
    hr = absTime(2);
    mn = absTime(3);
    sc = absTime(4);
    
    if((dowi==7)&&(dow==1))
        dd = dd+1;
    end
    DnTvec = [yyyy mm dd hr mn sc];
    DnTstr = datestr(DnTvec,'yyyy-mmm-dd, HH:MM:SS.FFF');
    fprintf(writeFileID, '%s, ', DnTstr);
    
    %leave rest of the columns of data unchanged
    for i = 5:13
        fprintf(writeFileID, '%s, ', cell2mat(readLine(i)));
    end
    fprintf(writeFileID, '%s\n', cell2mat(readLine(14)));
    
end

fclose(readFileID);
fclose(writeFileID);

%fprintf('attemtping to move file %s into %s \n', writeFile, readFile);
movefile(fullfile(dataDir,writeFile),fullfile(dataDir,readFile));

end





