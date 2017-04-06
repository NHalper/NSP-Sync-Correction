function CorrectDriftWithoutSignal (Method,DivideIntoThirds) 
%% Correct Drift Without Signal
% This function attempts to align two data files that have clock drift
% between the two. It attempts to determine the amount of clock drift by
% looking for cross correlations between various filtered versions of the
% signal between the two devices. If it fails to find a cross correlation
% that it feels confident in, it attempts to predict the clock drift based
% on emperical data for similar devices. 
% 
% Method:
% 0 = Script Chooses Best Method
% 1 = Force High Frequency Method
% 2 = Force Low Frequency Method
% 3 = Force Prediction Method
%
% DivideIntoThirds:
% 0 = Do not divide files into pieces
% 1 = Divide File Into Pieces
%
% Version 1.5
% - Corrected an issue with prediction method that would cause it %to not work with data that didn't contain cells.
%
% Version 1.6
% - Added a feature that allows one to select which NSP should be
% upsampled. 
% - Corrected an issue that would cause the program to fail to
% find NEV files under the low-frequency method.
%
% Version 1.7
% - Added modifications to process the files in pieces and therefore reduce
% RAM.
%
% Version 1.8
% - Removed the need for files to be processed in thirds.
% - Made the cross correlation measurements use a larger amount of data
% (300000 samples).
% - Made so that the lag calculations use averaged data.
% - Calculate lag at more points in the file and take a regression as the
% final lag calculation.
%
% Version 1.9
% - Allow division into thirds as an optional parameter
%
% Version 2.0
% - Allow use of paused files.
%


Version = '2.0';

TestMethod = 1; %Flag to set for faster processing method (basically now used as default method).
NumberOfFiles = 2; %Has to stay 2 for now as this script only operates on file pairs.
SkipThirds = 0; % A flag used later for small files. 
for i = 1:NumberOfFiles
    PacketLoss(i) = 0; % A flag set to attempt to correct packet loss. 
end

%% Select files to be drift corrected

uiwait(msgbox('Please select the files to align. One at a time.','Info','modal'));

% The NSP structure is repeatedly reused for containing either headers or
% actual full data structures. It isn't used in this block, but it is
% created here just to show its importance.
NSP = {};

for i = 1:NumberOfFiles   
    
    % This should no longer be required. Upsampled NSP is not selectable
    % anymore as the upsampled NSP is determined automatically in the
    % prediction method. (Upsampled NSP is otherwise determined by the lag
    % time in the xcorr calculation). Will remove this commented out line
    % in a future revision. 
%     if i == SelectUpsampleNSP
%         uiwait(msgbox('This NSP should be the 1.5. This one will be upsampled.','Info','modal'));
%     end
    
    % Open the file header without loading data. This lets us check the
    % sampling frequency and number (and identity) of channels in the file. 
    TempMetaData = openNSx('noread');
    
    %Require a 30k sampling for all files. 
    if TempMetaData.MetaTags.SamplingFreq ~= 30000
        disp('Signal is required to sampled at 30kHz. Please try again.')
        return
    end
    
    % Choose a random channel and use the next 5 channels. If there are
    % less than 5, then use all channels for the mean, but notify the user.
    % These channels are used to calculate the cross correlation to look
    % for high frequency events that are common to both. There is a channel
    % vector created per file loaded (NumberOfFiles)
    ChosenChannel = floor(rand*length(find([TempMetaData.ElectrodesInfo.ElectrodeID]<128)))+1;
    if length(find([TempMetaData.ElectrodesInfo.ElectrodeID]<128))>ChosenChannel+5
        Channel{i} = ChosenChannel:ChosenChannel+5;
    elseif length(find([TempMetaData.ElectrodesInfo.ElectrodeID]<128))>5
        Channel{i} = 1:5;
    else
        disp('Too few channels to take a 5 Channel mean signal. Using all channels.');
        Channel{i} = 1:length(find([TempMetaData.ElectrodesInfo.ElectrodeID]<128));
    end
    % Full filename used for openNSx later.
    Filename{i} = fullfile(TempMetaData.MetaTags.FilePath,[TempMetaData.MetaTags.Filename TempMetaData.MetaTags.FileExt]);
    clear ChosenChannel
end

%% Create a text based report of the results. Record version number. 

ReportFID = fopen(strcat(TempMetaData.MetaTags.FilePath,'/',TempMetaData.MetaTags.Filename,'-report.txt'),'w','n','Shift_JIS');
fprintf(ReportFID,'Report from RomanCorrection Script');
fprintf(ReportFID,'\n');
fprintf(ReportFID,'Version:');
fprintf(ReportFID,Version);
fprintf(ReportFID,'\n');

clear TempMetaData

%% Based on files selected above, open the channel subset for analysis

% Open NSx with channel parameter passed to open just five channels. Since
% they are consecutive, should only use five channels worth of memory for
% each entry in the NSP cell structure
for i = 1:NumberOfFiles
    NSP{i} = openNSx(Filename{i},['c:' num2str(Channel{i}(1)) ':' num2str(Channel{i}(end))]);
    NSPOrig{i} = openNSx(Filename{i},['c:' num2str(Channel{i}(1)) ':' num2str(Channel{i}(end))]);
end

% Save meta tags of original files into their own variable as the NSP
% variable is used up repeatedly. 
for i = 1:NumberOfFiles
    NSPMetaInfo{i} = openNSx(Filename{i},'noread');
end

clear Channel

%% Validate the data contents

% This whole loops checks the data for file recording restarts, resync
% events, and packet loss to make sure that it is operating on the right
% set of data. 
for i = 1:NumberOfFiles
    if iscell(NSP{i}.Data)
        disp('Cells in data');
        Cell = 1; % Flag to show that data is cellular
        SyncIndex(i) = FindReSync(NSP{i}.MetaTags); % External function finds clock resets. If it only finds one, it assumes that this is a resync event.
        if SyncIndex(i) == 0
            disp('File may have packet loss. This version of the script does not currently work with packet loss.')
            return
        end
        if SyncIndex(i) == -1 || SyncIndex(i) ~= length(NSP{i}.Data)
            disp('File may have packet loss. Will attempt to correct packet loss.')
            PacketLoss(i) = 1;
            if SyncIndex(i) == -1
                SyncIndex(i) = 1;
            end
        end
    else
        % Function can work without cell data, but script must change some
        % asssumptions to work this way. 
        disp('These files have single cell data.They may not have been from a synchronized recording, or they may have been processed by another tool.')
        Cell = 0;
    end
end

%% Test High Frequency Alignment Method
% This method attempts to look for a common high frequency signal between
% the files and use that signal to perform a cross correlation to determine
% the lag amount between the two files at that point. It validates this
% high frequency signal by looking for it in multiple places in the file
% and also comparing it to the known offset between the files.

% If no method is yet assigned or is the already assigned method, then this is considered a possible method and the data from this segment is gathered.  
if Method == 0 || Method == 1

    % Create the high pass filter used in this method.
    [b,a] = butter(8,2000/30000,'high');
    
    % Refresh the data
    NSP = NSPOrig;

    % Filter the data and remove portions that are not needed (pre-sync
    % periods)
    for i = 1:NumberOfFiles
        if Cell == 1
            if PacketLoss(i)
                for idx = 1:length(NSP{i}.Data)
                    NSP{i}.Data{idx} = filtfilt(b,a,double(NSP{i}.Data{idx})');
                end
            else  
                NSP{i}.Data = filtfilt(b,a,double(NSP{i}.Data{SyncIndex(i)})');
            end
        else
            NSP{i}.Data = filtfilt(b,a,double(NSP{i}.Data)');
        end
    end

    % Reduce dimensionality and noise by averaging the matrix we just filtered.
    for i = 1:NumberOfFiles
        if PacketLoss(i)
            for idx = 1:length(NSP{i}.Data)
                NSP{i}.Data{idx} = mean(NSP{i}.Data{idx}');
            end
        else
            if size(NSP{i}.Data,2) > 1
                NSP{i}.Data = mean(NSP{i}.Data');
            end
        end
    end

    % Find the shortest file time (can't compare files outside of that range)
    for i = 1:NumberOfFiles
        if PacketLoss(i)
            DataLength(i) = NSP{i}.MetaTags.Timestamp(end)+length(NSP{i}.Data{end});
        else
            DataLength(i) = length(NSP{i}.Data);
        end
    end
    EndPoint = min(DataLength);
    clear DataLength

    % Perform Cross Correlation
    for i = 1:NumberOfFiles
        if PacketLoss(i)
            DataToXcorr{i} = NSP{i}.Data{SyncIndex(i)}(1:300000);
        else
            DataToXcorr{i} = NSP{i}.Data(1:300000);
        end
    end
    
    SamplingRate = 30000; % Ensured 30k compliance earlier.
    [acor,lag] = xcorr(DataToXcorr{1}, DataToXcorr{2});
    [~, I] = max(abs(acor));
    firstlagdiff = lag(I); % Number of samples difference at this stage in the file. This is used for validation of our method by comparing it to timestamp offset.

    for i = 1:NumberOfFiles
        if Cell == 1
            if PacketLoss(i)
                DataTimestamps(i) = NSP{i}.MetaTags.Timestamp(SyncIndex(i));
            else
                DataTimestamps(i) = NSP{i}.MetaTags.Timestamp(SyncIndex(i)); %Get timestamp value of our data cell
            end
        else
            DataTimestamps(i) = NSP{i}.MetaTags.Timestamp; % If no cell, then only single timestamp value exists
        end
    end

    % Attempt to split file into thirds to save on RAM. These numbers
    % represent Channels to select when opening the files in the future. These
    % may not be used if the flag for passing the option to separate these into
    % thirds is not used. 
    for i = 1:NumberOfFiles
        NSPThirds{i,1} = floor(length(NSPMetaInfo{i}.MetaTags.ChannelID)/3);
        NSPThirds{i,2} = floor(2*(length(NSPMetaInfo{i}.MetaTags.ChannelID)/3));
        NSPThirds{i,3} = length(NSPMetaInfo{i}.MetaTags.ChannelID);
    end


    % If the high frequency method correctly predicts that actual offset of the
    % files then we have a winner. Allow some jitter/leniency.
    if abs(firstlagdiff) < (abs(max(DataTimestamps) - min(DataTimestamps))+10) || Method == 1
        disp('High frequency component found. Proceding with High Frequency method.')
        Method = 1; % 1==HighFrequencyMethod Later
        ModeComment = 'High Frequency';
    else
        disp('No high frequency component. Moving to test next method.')
        Method = 0; %Remain the same
    end

end


%% Test Low Frequency Alignment Method
% This method attempts to look for a common low frequency signal between
% the files and use that signal to perform a cross correlation to determine
% the lag amount between the two files at that point. It validates this
% low frequency signal by comparing it to the known offset between the
% files.

% If no method is yet assigned or is the already assigned method, then this is considered a possible method and the data from this segment is gathered.  
if Method == 0 || Method == 2

    % Create the high pass filter used in this method.
    [b,a] = butter(8,200/30000,'low');
    
    % Refresh the data
    NSP = NSPOrig;

    % Filter the data and remove portions that are not needed (pre-sync
    % periods)
    for i = 1:NumberOfFiles
        if Cell == 1
            if PacketLoss(i)
                for idx = 1:length(NSP{i}.Data)
                    NSP{i}.Data{idx} = filtfilt(b,a,double(NSP{i}.Data{idx})');
                end
            else  
                NSP{i}.Data = filtfilt(b,a,double(NSP{i}.Data{SyncIndex(i)})');
            end
        else
            NSP{i}.Data = filtfilt(b,a,double(NSP{i}.Data)');
        end
    end

    % Reduce dimensionality and noise by averaging the matrix we just filtered.
    for i = 1:NumberOfFiles
        if PacketLoss(i)
            for idx = 1:length(NSP{i}.Data)
                NSP{i}.Data{idx} = mean(NSP{i}.Data{idx}');
            end
        else
            if size(NSP{i}.Data,2) > 1
                NSP{i}.Data = mean(NSP{i}.Data');
            end
        end
    end

    % Find the shortest file time (can't compare files outside of that range)
    for i = 1:NumberOfFiles
        if PacketLoss(i)
            DataLength(i) = NSP{i}.MetaTags.Timestamp(end)+length(NSP{i}.Data{end});
        else
            DataLength(i) = length(NSP{i}.Data);
        end
    end
    EndPoint = min(DataLength);
    clear DataLength

    % Perform Cross Correlation
    for i = 1:NumberOfFiles
        if PacketLoss(i)
            DataToXcorr{i} = NSP{i}.Data{SyncIndex(i)}(1:300000);
        else
            DataToXcorr{i} = NSP{i}.Data(1:300000);
        end
    end
    
    SamplingRate = 30000; % Ensured 30k compliance earlier.
    [acor,lag] = xcorr(DataToXcorr{1}, DataToXcorr{2});
    [~, I] = max(abs(acor));
    firstlagdiff = lag(I); % Number of samples difference at this stage in the file. This is used for validation of our method by comparing it to timestamp offset.

    for i = 1:NumberOfFiles
        if Cell == 1
            if PacketLoss(i)
                DataTimestamps(i) = NSP{i}.MetaTags.Timestamp(SyncIndex(i));
            else
                DataTimestamps(i) = NSP{i}.MetaTags.Timestamp(SyncIndex(i)); %Get timestamp value of our data cell
            end
        else
            DataTimestamps(i) = NSP{i}.MetaTags.Timestamp; % If no cell, then only single timestamp value exists
        end
    end

    % Attempt to split file into thirds to save on RAM. These numbers
    % represent Channels to select when opening the files in the future. These
    % may not be used if the flag for passing the option to separate these into
    % thirds is not used. 
    for i = 1:NumberOfFiles
        NSPThirds{i,1} = floor(length(NSPMetaInfo{i}.MetaTags.ChannelID)/3);
        NSPThirds{i,2} = floor(2*(length(NSPMetaInfo{i}.MetaTags.ChannelID)/3));
        NSPThirds{i,3} = length(NSPMetaInfo{i}.MetaTags.ChannelID);
    end


    % If the low frequency method correctly predicts that actual offset of the
    % files then we have a winner. Allow some jitter/leniency.
    if abs(firstlagdiff) < (abs(max(DataTimestamps) - min(DataTimestamps))+10) || Method == 2
        disp('Low frequency component found. Proceding with Low Frequency method.')
        Method = 2; % 2==LowFrequencyMethod
        ModeComment = 'Low Frequency';
    else
        disp('No low frequency component. Moving to test next method.')
        Method = 0; %Remain the same
    end
end


%% Test Prediction Alignment Method
% This method uses empirical data for similar devices to determine what the
% amount of offset should be over a given time period. 

% If no method is yet assigned or is the already assigned method, then this is considered a possible method and the data from this segment is gathered.  
if Method == 0 || Method == 3
    
    % Refresh the Data
    NSP = NSPOrig;
    
    % Find the shortest file time (can't compare files outside of that range)
    for i = 1:NumberOfFiles
        if PacketLoss(i)
            DataLength(i) = NSP{i}.MetaTags.Timestamp(end)+length(NSP{i}.Data{end})-NSP{i}.MetaTags.Timestamp(SyncIndex(i));
        else
            DataLength(i) = length(NSP{i}.Data);
        end
    end
    EndPoint = min(DataLength);
    clear DataLength
    
    % This will be used to determine the number of analog inputs in each
    % file. This is used to determine which device is the one that needs to
    % be upsampled according to empircal data.
    AnalogInputCounts = [];
    for i = 1:NumberOfFiles
        AnalogInputCounts(i) = length(find([NSP{i}.ElectrodesInfo.ElectrodeID]>128));
    end

    %Check for event codes in digital inputs of NEV file. This will replace
    %the section below
    DataToResample = [];
     for i = 1:NumberOfFiles
         if (exist(fullfile(NSP{i}.MetaTags.FilePath, [NSP{i}.MetaTags.Filename '.nev'])))
             TempNEV = openNEV(fullfile(NSP{i}.MetaTags.FilePath, [NSP{i}.MetaTags.Filename '.nev']),'nosave','nomat');
             disp(strcat('Found:',fullfile(NSP{i}.MetaTags.FilePath, [NSP{i}.MetaTags.Filename '.nev'])));
             if ~isempty(TempNEV.Data.SerialDigitalIO.TimeStamp) 
                 if isempty(DataToResample)
                     DataToResample = i;
                 else
                     disp('Multiple files contain digital events. Cannot distinguish NSP 1.5');
                     return
                 end
             end
         end 
     end
     
     if isempty(DataToResample)
         disp('No digital events found (or no NEV files found), so script cannot determine upsample NSP using prediction method.')
         disp('Hard code which NSP should be chosen by editing this value below, or make sure that the NEV files are visible to the script.')
         DataToResample = 1;
         %return
     end
     
    % Filter the data and remove portions that are not needed (pre-sync
    % periods)
    for i = 1:NumberOfFiles
        if Cell == 1
            if PacketLoss(i) % Remove pauses.
                NumChan = size(NSP{i}.Data{1},1);
                for ipq = SyncIndex(i):length(NSP{i}.Data)-1
                    PauseTimestamp = NSP{i}.MetaTags.Timestamp(ipq) + length(NSP{i}.Data{ipq});
                    NextSegTimestamp = NSP{i}.MetaTags.Timestamp(ipq+1);
                    NSP{i}.Data{ipq} = [NSP{i}.Data{ipq} zeros(NumChan,NextSegTimestamp-PauseTimestamp)];
                end
                NSP{i}.Data = cell2mat(NSP{i}.Data(SyncIndex(i):end))';
                NSP{i}.MetaTags.Timestamp = NSP{i}.MetaTags.Timestamp(SyncIndex(i));
            else  
                NSP{i}.Data = NSP{i}.Data{SyncIndex(i)}';
            end
        else
            NSP{i}.Data = NSP{i}.Data;
        end
    end

    % Reduce dimensionality and noise by averaging the matrix we just filtered.
    for i = 1:NumberOfFiles
        if size(NSP{i}.Data,1) > 1
            NSP{i}.Data = mean(NSP{i}.Data');
        end
    end

    % Perform cross correlation
    SamplingRate = 30000; % Ensured 30k compliance earlier.
    [acor,lag] = xcorr(NSP{1}.Data(1:300000), NSP{2}.Data(1:300000));
    [~, I] = max(abs(acor));
    firstlagdiff = lag(I); % Number of samples difference at this stage in the file. This is used for validation of our method by comparing it to timestamp offset.

    for i = 1:NumberOfFiles
        DataTimestamps(i) = NSP{i}.MetaTags.Timestamp; % If no cell, then only single timestamp value exists
    end

    % Attempt to split file into thirds to save on RAM. These numbers
    % represent Channels to select when opening the files in the future. These
    % may not be used if the flag for passing the option to separate these into
    % thirds is not used. 
    for i = 1:NumberOfFiles
        NSPThirds{i,1} = floor(length(NSPMetaInfo{i}.MetaTags.ChannelID)/3);
        NSPThirds{i,2} = floor(2*(length(NSPMetaInfo{i}.MetaTags.ChannelID)/3));
        NSPThirds{i,3} = length(NSPMetaInfo{i}.MetaTags.ChannelID);
    end


    % Since this is the last option. Force prediction method.
    disp('No frequency component found. Proceding with Prediction method.')
    Method = 3; % 3==LowFrequencyMethod Later
    ModeComment = 'Prediction Method';

end


%% Apply the Given Method

% Add method information to the report
fprintf(ReportFID,'Mode:');
fprintf(ReportFID,ModeComment);
fprintf(ReportFID,'\n');

% Perform Cross Correlation
for i = 1:NumberOfFiles
    if PacketLoss(i)
        if Method == 3
            DataToXcorr{i} = NSP{i}.Data(EndPoint-300000:EndPoint);
        else
            DataToXcorr{i} = NSP{i}.Data{end}(EndPoint-NSP{i}.MetaTags.Timestamp(end)-300000:EndPoint-NSP{i}.MetaTags.Timestamp(end));
        end
    else
        DataToXcorr{i} = NSP{i}.Data(EndPoint-300000:EndPoint);
    end
end

if PacketLoss(1) == 1 || PacketLoss(2) == 1
    [acor,lag] = xcorr(DataToXcorr{1}, DataToXcorr{2});
    [~, I] = max(abs(acor));
    EndLag = lag(I);
else
    % Perform Cross Correlation at various points in the file as a way of
    % validating the signal 
    lagdiff = [];% Used to store the lag of the file near the endpoint. Difference between this and first lagdiff is total drift up to that point in the file.
    NumberOfLagPoints = 20;
    for i = 1:NumberOfLagPoints
        if i == NumberOfLagPoints
            [acor,lag] = xcorr(NSP{1}.Data(EndPoint-300000:EndPoint), NSP{2}.Data(EndPoint-300000:EndPoint));
        else
            LagSamplePoint = round(EndPoint/(NumberOfLagPoints-i));
            if LagSamplePoint < 300000
                [acor,lag] = xcorr(NSP{1}.Data(1:LagSamplePoint), NSP{2}.Data(1:LagSamplePoint));
            else
                [acor,lag] = xcorr(NSP{1}.Data(LagSamplePoint-300000:LagSamplePoint), NSP{2}.Data(LagSamplePoint-300000:LagSamplePoint));
            end
        end
        [~, I] = max(abs(acor));
        %plot(acor);
        lagdiff(i) = lag(I);
    end
    clear LagSamplePoint
    clear SamplingRate
    clear acor
    clear lag
    clear I

    lm = fitlm([1:length(lagdiff)],lagdiff);
    Intercept = lm.Coefficients{1,1};
    Coefficient = lm.Coefficients{2,1};
    EndLag = round(Intercept+Coefficient*NumberOfLagPoints);
    LagDifferenceCalculation = abs(lagdiff(NumberOfLagPoints)-EndLag);

    fprintf(ReportFID,'Intercept: ');
    fprintf(ReportFID,num2str(Intercept));
    fprintf(ReportFID,'\n');

    fprintf(ReportFID,'Coefficient: ');
    fprintf(ReportFID,num2str(Coefficient));
    fprintf(ReportFID,'\n');

    fprintf(ReportFID,'R Squared: ');
    fprintf(ReportFID,num2str(lm.Rsquared.Ordinary));
    fprintf(ReportFID,'\n');
end

% If Lag Amount is Negative, NSP1 Data needs to be resampled at a higher
% rate. If Positive, NSP2 needs to be resampled at a higher rate.
if Method == 3
    %Data to Resample is already set.
else
    if (lagdiff(1) < 0)
        DataToResample = 1;
    elseif (lagdiff(1) > 0)
        DataToResample = 2;
    else
        disp('Your drift amount appears to be 0; these data may be aligned. Please check manually')
        return
    end
end

% Original length is/was used as a check to see that our total data
% length had increased by the expected amount as a crude way of
% checking if resampling had been properly aplied. 
if PacketLoss(DataToResample)
    if Method == 3
        OriginalLength = length(NSP{DataToResample}.Data);
    else
        OriginalLength = length(NSP{DataToResample}.Data{end}) + NSP{DataToResample}.MetaTags.Timestamp(end);
    end
else
    OriginalLength = length(NSP{DataToResample}.Data);
end

% Add to report the timing offset throughout the file.
if PacketLoss(1) == 1 || PacketLoss(2) == 1
    fprintf(ReportFID,['lagdiff: ']);
    fprintf(ReportFID,num2str(EndLag));
    fprintf(ReportFID,'\n');
else
    for i = 1:NumberOfLagPoints
        fprintf(ReportFID,['lagdiff ' num2str(i) ': ']);
        fprintf(ReportFID,num2str(lagdiff(i)));
        fprintf(ReportFID,'\n');
    end
end
% Total drift amount is one of the most useful values generated. It
% is saved to the report here. 
if Method == 3
    DriftAmount = 100*(OriginalLength/(30000*60*5)); %We see around 100 samples loss every 5 minutes on most of these NSPs.
else
    DriftAmount = abs(EndLag - firstlagdiff);
end
fprintf(ReportFID,'Total Drift Amount:');
fprintf(ReportFID,num2str(DriftAmount));
fprintf(ReportFID,'\n');
clear NumberOfLagPoints

% Resample period tells us every how many samples a new sample
% needs to be added. 
ResamplePeriod = round(EndPoint/DriftAmount);

% Begin to work on the data. 
disp('Calculations complete. Opening full data file for drift correction. This may take a long while.')

% Clear NSP structure since it will be redefined in the next functional
% line. This is less time efficient, but seems to be more memory
% efficient. It may be bad practice to reuse a variable in a very
% different way, but since the NSP structure is a way to maintain
% consistency in information type in this script, I did it.
clear NSP
clear NSPOrig

% Do once for each file to resample the original 30k data. This will
% involve opening the whole file so could be very RAM intensive. 
for i = 1:NumberOfFiles
    for idx = 1:(1 + DivideIntoThirds*2)
        if DivideIntoThirds == 1
            switch(idx)
                case 1
                    NSP{i} = openNSx(Filename{i},['c:' num2str(1) ':' num2str(NSPThirds{i,1})]);
                case 2
                    NSP{i} = openNSx(Filename{i},['c:' num2str(NSPThirds{i,1}+1) ':' num2str(NSPThirds{i,2})]);
                case 3
                    NSP{i} = openNSx(Filename{i},['c:' num2str(NSPThirds{i,2}+1) ':' num2str(NSPThirds{i,3})]);
            end
        else
            NSP{i} = openNSx(Filename{i});
        end

        % If data is cellular (has proper resync events), only a portion of it is needed. 
        if Cell == 1     
            if PacketLoss(i)
                NumChan = size(NSP{i}.Data{1},1);
                for ipq = SyncIndex(i):length(NSP{i}.Data)-1
                    PauseTimestamp = NSP{i}.MetaTags.Timestamp(ipq) + length(NSP{i}.Data{ipq});
                    NextSegTimestamp = NSP{i}.MetaTags.Timestamp(ipq+1);
                    NSP{i}.Data{ipq} = [NSP{i}.Data{ipq} zeros(NumChan,NextSegTimestamp-PauseTimestamp)];
                end
                NSP{i}.Data = cell2mat(NSP{i}.Data(SyncIndex(i):end));
                NSP{i}.MetaTags.Timestamp = NSP{i}.MetaTags.Timestamp(SyncIndex(i));
            else  
                NSP{i}.Data = NSP{i}.Data{SyncIndex(i)};
                NSP{i}.MetaTags.Timestamp = NSP{i}.MetaTags.Timestamp(SyncIndex(i));
            end
        else
            % No cells in data, so it can be used as is. 
        end    

        TotalPeriods = round(EndPoint/ResamplePeriod); %This is equal to drift amount, but it makes the logic a little bit easier to follow for new people.

        % This is where the actual resampling happens.
        % This is the fastest and most efficient method I could
        % find for repeating given elements of the array as a
        % method of upsampling. 
        if i == DataToResample
            tic;
            RepeatingArray = ones(1,length(NSP{DataToResample}.Data));
            RepeatingArray(1:ResamplePeriod:end) = 2;
            NSP{DataToResample}.Data = repelem(NSP{DataToResample}.Data,1,RepeatingArray);
            toc;
        end

        % Correct disparity between the files caused by resynchronization.
        % If we did this earlier, then the xcorr could have accidentally
        % latched on to this manual offset, doing this now just reduces the
        % risk of that happening and ensures we have the expected offset
        % when validating this method. 
        NSP{i}.Data = [zeros(size(NSP{i}.Data,1),NSP{i}.MetaTags.Timestamp) NSP{i}.Data];
        NSP{i}.MetaTags.Timestamp = 0;

        % Save the files using the saveNSxSync function, which is like
        % SaveNSx, but doesn't have user input prompts. 
        saveNSxSync(NSP{i},0);

        % It is about to be redefined again in the next iteration of the
        % loop. If it is the last iteration, then we don't need it anymore!
        clear NSP
    end
end

clear NSPThirds

%% Additional File Types
% Handle Additional Continous File Types (This is how the loop
% structure should be for above, but above the nesting is reversed for
% some reason. Fixing that would allow one to toggle between splitting the file into thirds or
% not...)
for i = 1:NumberOfFiles

    % Creates a matrix of file types and elminates ones that have
    % already been used. 
    if PacketLoss(i) ~= 1
        ValidFileTypes = {};
            for idx = 1:NumberOfFiles
                FileExtTypes = {'.ns1' '.ns2' '.ns3' '.ns4' '.ns5' '.ns6'};
                ValidTypesIndices = [];
                for Type = 1:length(FileExtTypes)
                    if ~strcmpi(FileExtTypes{Type},NSPMetaInfo{i}.MetaTags.FileExt)
                        ValidTypesIndices = [ValidTypesIndices Type];
                    end
                end
                ValidFileTypes{idx} = ValidTypesIndices;
            end

        % For each available file type, it attempts to look for a file.
        for Type = ValidFileTypes{i}
            TempFilename = fullfile(NSPMetaInfo{i}.MetaTags.FilePath, [NSPMetaInfo{i}.MetaTags.Filename FileExtTypes{Type}]);

            % If the file exists, it is time to perform the exact same
            % things that we did above, but on the new file. 
            if exist(TempFilename)
                disp(strcat('Found:',TempFilename));
                TempStructure = openNSx(TempFilename,'noread');
                EndPacket = TempStructure.MetaTags.DataPoints;
                TimestampScale = 30000/TempStructure.MetaTags.SamplingFreq; % Used to scale how many timestamps to insert (downsampled upsample value)
                if DivideIntoThirds == 1
                    SubSectionRepeats = 3;
                else
                    SubSectionRepeats = 1;
                end
                for ixxy = 1:SubSectionRepeats
                    NSPThirds(1) = floor(length(TempStructure.MetaTags.ChannelID)/3);
                    if NSPThirds(1) == 0
                        NSPThirds(1) = 1;
                    end
                    NSPThirds(2) = floor(2*(length(TempStructure.MetaTags.ChannelID)/3));
                    NSPThirds(3) = length(TempStructure.MetaTags.ChannelID);

                    if range(NSPThirds) < 2
                        SkipThirds = 1;
                    end

                    if DivideIntoThirds == 1 && SkipThirds == 0 % If divide into thirds, then ixxy will be 1:3 and this switch will hit each third, otherwise skip and just open the file. ixxy will be 1 and this loop will go once.
                        switch(ixxy)
                            case 1;
                                NSx = openNSx(TempFilename,['c:' num2str(1) ':' num2str(NSPThirds(1))]);
                            case 2;
                                NSx = openNSx(TempFilename,['c:' num2str(NSPThirds(1)+1) ':' num2str(NSPThirds(2))]); 
                            case 3;
                                NSx = openNSx(TempFilename,['c:' num2str(NSPThirds(2)+1) ':' num2str(NSPThirds(3))]); 
                        end
                    else
                        NSx = openNSx(TempFilename); 
                    end

                    if DataToResample == i
                        if and(and(iscell(NSx.Data),Cell==1),length(NSx.Data)==2)
                            NSx.Data = [zeros(size(NSx.Data{2},1),round(NSx.MetaTags.Timestamp(2)/30000)) NSx.Data{2}];
                            NSx.MetaTags.Timestamp = 0;
                            TempResamplePeriod = round((ResamplePeriod/TimestampScale)*TimestampScale);
                            TempEndPoint = round(EndPoint/TimestampScale);
                            TotalPeriods = floor(TempEndPoint/TempResamplePeriod);
                            if TestMethod == 1
                                tic
                                RepeatingArray = ones(1,length(NSx.Data));
                                RepeatingArray(1:TempResamplePeriod:end) = 2;
                                NSx.Data = repelem(NSx.Data,1,RepeatingArray);
                                toc
                            else
                                if TotalPeriods > 1
                                    for period = 1:TotalPeriods
                                        NSx.Data = [NSx.Data(:,1:TempResamplePeriod*period) NSx.Data(:,TempResamplePeriod*period:end)];
                                        disp([num2str(period) 'of' num2str(TotalPeriods)]);
                                    end
                                end
                            end

                        elseif not(iscell(NSx.Data))
                            NSx.Data = [zeros(size(NSx.Data,1),NSx.MetaTags.Timestamp) NSx.Data];
                            NSx.MetaTags.Timestamp = 0;
                            TempResamplePeriod = round((ResamplePeriod/TimestampScale)*TimestampScale);
                            TempEndPoint = round(EndPoint/TimestampScale);
                            TotalPeriods = floor(TempEndPoint/TempResamplePeriod);
                            if TestMethod == 1
                                tic
                                RepeatingArray = ones(1,length(NSx.Data));
                                RepeatingArray(1:TempResamplePeriod:end) = 2;
                                NSx.Data = repelem(NSx.Data,1,RepeatingArray);
                                toc
                            else
                                if TotalPeriods > 1
                                    for period = 1:TotalPeriods
                                        NSx.Data = [NSx.Data(:,1:TempResamplePeriod*period) NSx.Data(:,TempResamplePeriod*period:end)];
                                        disp([num2str(period) 'of' num2str(TotalPeriods)]);
                                    end
                                end
                            end
                        end
                    end
                    saveNSxSync(NSx,0);
                    if range(NSPThirds) < 2
                        break;
                    end
                end
            end
        end
    end
    % We are done with continuous data segments, so we can clear this
    % variable to save on some memory.
    clear NSx 

%% Handle NEV Resampling

    % This is still part of the 'for number of files' loop above, but
    % it handles NEV instead of continuous data. The base method for
    % this is to properly add values to each timestamp based on how far
    % into the file it is. This means that we must iterate over each
    % timestamp since each one will have a potentially different value
    % added to it depending on its own value. Overall, though, the
    % looping is pretty fast even for long files as they are 1D arrays
    % and we are just doing addition. 
    if (exist(fullfile(NSPMetaInfo{i}.MetaTags.FilePath, [NSPMetaInfo{i}.MetaTags.Filename '.nev'])))
        TempNEV = openNEV(fullfile(NSPMetaInfo{i}.MetaTags.FilePath, [NSPMetaInfo{i}.MetaTags.Filename '.nev']),'nosave','nomat');
        disp(strcat('Found:',fullfile(NSPMetaInfo{i}.MetaTags.FilePath, [NSPMetaInfo{i}.MetaTags.Filename '.nev'])));
        %Serial Digital IO
        for timestamp = 1:length(TempNEV.Data.SerialDigitalIO.TimeStamp)
            TempNEV.Data.SerialDigitalIO.TimeStamp(timestamp) = TempNEV.Data.SerialDigitalIO.TimeStamp(timestamp)+round(TempNEV.Data.SerialDigitalIO.TimeStamp(timestamp)/ResamplePeriod);
        end
        %Spikes
        for timestamp = 1:length(TempNEV.Data.Spikes.TimeStamp)
            TempNEV.Data.Spikes.TimeStamp(timestamp) = TempNEV.Data.Spikes.TimeStamp(timestamp)+round(TempNEV.Data.Spikes.TimeStamp(timestamp)/ResamplePeriod);
        end
        %Comments
        for timestamp = 1:length(TempNEV.Data.Comments.TimeStamp)
            TempNEV.Data.Comments.TimeStamp(timestamp) = TempNEV.Data.Comments.TimeStamp(timestamp)+round(TempNEV.Data.Comments.TimeStamp(timestamp)/ResamplePeriod);
        end
        %Video Sync
        for timestamp = 1:length(TempNEV.Data.VideoSync.TimeStamp)
            TempNEV.Data.VideoSync.TimeStamp(timestamp) = TempNEV.Data.VideoSync.TimeStamp(timestamp)+round(TempNEV.Data.VideoSync.TimeStamp(timestamp)/ResamplePeriod);
        end
        %Tracking Events
        for timestamp = 1:length(TempNEV.Data.TrackingEvents.TimeStamp)
            TempNEV.Data.TrackingEvents.TimeStamp(timestamp) = TempNEV.Data.TrackingEvents.TimeStamp(timestamp)+round(TempNEV.Data.TrackingEvents.TimeStamp(timestamp)/ResamplePeriod);
        end
        %Patient Trigger
        for timestamp = 1:length(TempNEV.Data.PatientTrigger.TimeStamp)
            TempNEV.Data.PatientTrigger.TimeStamp(timestamp) = TempNEV.Data.PatientTrigger.TimeStamp(timestamp)+round(TempNEV.Data.PatientTrigger.TimeStamp(timestamp)/ResamplePeriod);
        end
        %Reconfig
        for timestamp = 1:length(TempNEV.Data.Reconfig.TimeStamp)
            TempNEV.Data.Reconfig.TimeStamp(timestamp) = TempNEV.Data.Reconfig.TimeStamp(timestamp)+round(TempNEV.Data.Reconfig.TimeStamp(timestamp)/ResamplePeriod);
        end

        if ~isempty(TempNEV.Data.Tracking)
            uiwait(msgbox('This version does not operate on tracking data timestamps. Please contact Blackrock Support.','ERROR','modal'));
        end

        saveNEVSync(TempNEV);
    end
end

% Finish the report.
fclose(ReportFID);

%% Wrap Up
return

