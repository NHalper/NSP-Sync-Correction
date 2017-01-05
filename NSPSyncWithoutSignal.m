function NSPSyncWithoutSignal (ForceHFMethod,ForceLFMethod,ForcePredictionMethod,SelectUpsampleNSP) 
%% 
% 
% 
% 
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
Version = '1.6';

TestMethod = 1; %Flag to set for faster processing method

NumberOfFiles = 2; %Has to stay 2 for now

Method = 0; %This is a flag that will be set when a certain method of alignment should be used.

uiwait(msgbox('Please select the files to align. One at a time.','Info','modal'));

NSP = {};
for i = 1:NumberOfFiles   
    if i == SelectUpsampleNSP
        uiwait(msgbox('This NSP should be the 1.5. This one will be upsampled.','Info','modal'));
    end
    TempMetaData = openNSx('noread');
    if TempMetaData.MetaTags.SamplingFreq ~= 30000
        disp('Signal is required to sampled at 30kHz. Please try again.')
        return
    end
    %Choosing a random channel
    Channel{i} = num2str(floor(rand*length(find([TempMetaData.ElectrodesInfo.ElectrodeID]<128)))+1);
    Filename{i} = fullfile(TempMetaData.MetaTags.FilePath,[TempMetaData.MetaTags.Filename TempMetaData.MetaTags.FileExt]);
end

ReportFID = fopen(strcat(TempMetaData.MetaTags.FilePath,'/',TempMetaData.MetaTags.Filename,'-report.txt'),'w','n','Shift_JIS');
fprintf(ReportFID,'Report from RomanCorrection Script');
fprintf(ReportFID,'\n');
fprintf(ReportFID,'Version:');
fprintf(ReportFID,Version);
fprintf(ReportFID,'\n');

for i = 1:NumberOfFiles
    NSP{i} = openNSx(Filename{i},strcat('c:',Channel{i},':',Channel{i}));
end
% clear TempMetaData
clear Channel

for i = 1:NumberOfFiles
    if iscell(NSP{i}.Data)
        disp('Cells in data');
        Cell = 1;
        SyncIndex(i) = FindReSync(NSP{i}.MetaTags);
        if SyncIndex(i) == 0
            return
        end
        if SyncIndex(i) < length(NSP{i}.Data)-1
            disp('File may have packet loss. This version of the script does not currently work with packet loss.')
            return
        end
    else
        disp('These files have single cell data.They may not have been from a synchronized recording, or they may have been processed by another tool.')
        Cell = 0;
    end
end


%% Attempt to Look for a common High Frequency Signal

%Create the high pass filter.
[b,a] = butter(8,2000/30000,'high');

%Filter the data and remove portions that are not needed
for i = 1:NumberOfFiles
    if Cell == 1
        NSP{i}.Data = filtfilt(b,a,double(NSP{i}.Data{SyncIndex(i)})');
    else
        NSP{i}.Data = filtfilt(b,a,double(NSP{i}.Data)');
    end
end

%Correct disparity between the files caused by resynchronization
% for i = 1:NumberOfFiles
%     NSP{i}.Data = [zeros(size(NSP{i}.Data,1),NSP{i}.MetaTags.Timestamp(SyncIndex)) NSP{i}.Data];
% end

%Perform Cross Correlation
lagdiff = [];
SamplingRate = 30000;
[acor,lag] = xcorr(NSP{1}.Data(1:30000), NSP{2}.Data(1:30000));
[~, I] = max(abs(acor));
firstlagdiff = lag(I)
% timediff = lagdiff/Fs

for i = 1:NumberOfFiles
    if Cell == 1;
        DataTimestamps(i) = NSP{i}.MetaTags.Timestamp(SyncIndex(i));
    else
        DataTimestamps(i) = NSP{i}.MetaTags.Timestamp;
    end
end

% If the high frequency method correctly predicts that actual offset of the
% files then we have a winner. Allow some jitter/leniency
if abs(firstlagdiff) < (abs(max(DataTimestamps) - min(DataTimestamps))+10)
    disp('High frequency component found. Proceding with High Frequency method.')
    Method = 1;
else
    disp('No high frequency component. Moving to next method.')
    Method = 0; %Remain the same
end

if ForceHFMethod == 1
    Method = 1;
end

if ForceLFMethod == 1
    Method = 0;
end

if ForcePredictionMethod == 1
    Method = 0;
end


if Method == 1
    fprintf(ReportFID,'Mode:');
    fprintf(ReportFID,'High Frequency');
    fprintf(ReportFID,'\n');
    
    for i = 1:NumberOfFiles
        DataLength(i) = length(NSP{i}.Data);
    end
    EndPoint = min(DataLength);
    clear DataLength

    %Perform Cross Correlation
    lagdiff = [];
    SamplingRate = 30000;
    [acor,lag] = xcorr(NSP{1}.Data(EndPoint-30000:EndPoint), NSP{2}.Data(EndPoint-30000:EndPoint));
    [~, I] = max(abs(acor));
    plot(acor);
    lagdiff = lag(I)
    clear SamplingRate
    clear acor
    clear lag
    clear I


    %If Lag Amount is Negative, NSP1 Data needs to be resampled at a higher
    %rate. If Positive, NSP2 needs to be resampled at a higher rate.
    if (lagdiff < 0)
        DataToResample = 1;
    elseif (lagdiff > 0)
        DataToResample = 2;
    else
        disp('Your drift amount appears to be 0; these data may be aligned. Please check manually')
        DataToResample = 0;
        return
    end

    disp('Calculations complete. Opening full data file for drift correction. This may take a long while.')
    fprintf(ReportFID,'lagdiff:');
    fprintf(ReportFID,num2str(lagdiff));
    fprintf(ReportFID,'\n');
    for i = 1:NumberOfFiles
        NSP{i} = openNSx(Filename{i});
        if Cell == 1;     
            NSP{i}.Data = NSP{i}.Data{SyncIndex(i)};
            NSP{i}.MetaTags.Timestamp = NSP{i}.MetaTags.Timestamp(SyncIndex(i));
        else

        end
    end
    DriftAmount = abs(lagdiff - firstlagdiff)
    fprintf(ReportFID,'Total Drift Amount:');
    fprintf(ReportFID,num2str(DriftAmount));
    fprintf(ReportFID,'\n');
    OriginalLength = length(NSP{DataToResample}.Data);
    ResamplePeriod = round(EndPoint/DriftAmount);
%     OriginalSamplingRate = 30000;
%     NewSamplingRate = 30000 + 30000/ResamplePeriod;
%     [p,q] = rat(NewSamplingRate/OriginalSamplingRate);
%     DataForResampling = resample(double(DataForResampling),p,q);
    TotalPeriods = round(EndPoint/ResamplePeriod); %This is basically just equal to drift amount, but it looks nice
    
    if TestMethod == 1
        tic
            RepeatingArray = ones(1,length(NSP{DataToResample}.Data));
            RepeatingArray(1:ResamplePeriod:end) = 2;
            NSP{DataToResample}.Data = repelem(NSP{DataToResample}.Data,1,RepeatingArray);
        toc
    else
        for period = 1:TotalPeriods
                tic
                NSP{DataToResample}.Data = [NSP{DataToResample}.Data(:,1:ResamplePeriod*period) NSP{DataToResample}.Data(:,ResamplePeriod*period:end)];
                disp([num2str(period) 'of' num2str(TotalPeriods)]);
                t = toc;
                EstimatedTimeLeft = t*(TotalPeriods-period);
                disp(strcat('Estimated Time Left: ',num2str(EstimatedTimeLeft),' Seconds'))
        end
    end
    
    %Correct disparity between the files caused by resynchronization
    for i = 1:NumberOfFiles
        NSP{i}.Data = [zeros(size(NSP{i}.Data,1),NSP{i}.MetaTags.Timestamp) NSP{i}.Data];
    end
    
    MethodComment = 'Sync Method: High Frequency Method';
    NSP{1}.MetaTags.Comment(end-length(MethodComment)+1:end) = MethodComment;
    NSP{2}.MetaTags.Comment(end-length(MethodComment)+1:end) = MethodComment;
    saveNSxSync(NSP{1});
    saveNSxSync(NSP{2});
    
    %Handle Additional File Types
    for i = 1:NumberOfFiles
        ValidFileTypes = {};
            for idx = 1:NumberOfFiles
                FileExtTypes = {'.ns1' '.ns2' '.ns3' '.ns4' '.ns5' '.ns6'};
                ValidTypesIndices = [];
                for Type = 1:length(FileExtTypes)
                    if ~strcmpi(FileExtTypes{Type},NSP{1}.MetaTags.FileExt)
                        ValidTypesIndices = [ValidTypesIndices Type];
                    end
                end
                ValidFileTypes{idx} = ValidTypesIndices;
            end
            
            
        for Type = ValidFileTypes{i}
            TempFilename = fullfile(NSP{i}.MetaTags.FilePath, [NSP{i}.MetaTags.Filename FileExtTypes{Type}])
            
            if exist(TempFilename)
                disp(strcat('Found:',TempFilename));
                TempStructure = openNSx(TempFilename,'noread');
                EndPacket = TempStructure.MetaTags.DataPoints;
                TimestampScale = 30000/TempStructure.MetaTags.SamplingFreq;
                clear TempStructure
                NSx = openNSx(TempFilename);
                if DataToResample == i
                    if and(and(iscell(NSx.Data),Cell==1),length(NSx.Data)==2)
                        %disp('This function does not work on paused data currently');

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
                        TempResamplePeriod = round((ResamplePeriod/TimestampScale)*TimestampScale);
                        TempEndPoint = round(EndPoint/TimestampScale);
                        TotalPeriods = floor(TempEndPoint/TempResamplePeriod);
                        if TotalPeriods > 1
                            for period = 1:TotalPeriods
                                NSx.Data = [NSx.Data(:,1:TempResamplePeriod*period) NSx.Data(:,TempResamplePeriod*period:end)];
                                disp([num2str(period) 'of' num2str(TotalPeriods)]);
                            end
                        end
                    end
                end
                
                
                saveNSxSync(NSx);
            end
        end
        
        
        %NEV Attempt
        if (exist(fullfile(NSP{i}.MetaTags.FilePath, [NSP{i}.MetaTags.Filename '.nev'])))
            TempNEV = openNEV(fullfile(NSP{i}.MetaTags.FilePath, [NSP{i}.MetaTags.Filename '.nev']),'nosave','nomat');
            disp(strcat('Found:',fullfile(NSP{i}.MetaTags.FilePath, [NSP{i}.MetaTags.Filename '.nev'])));
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
    
    
    fclose(ReportFID);
    return
else
    %disp('No high frequency component moving on to low frequency method.')
end 

%% Attempt to look for common low frequency signal

%Create the low pass filter.
[b,a] = butter(8,200/30000,'low');

%Filter the data and remove portions that are not needed
for i = 1:NumberOfFiles
    if Cell == 1
        NSP{i}.Data = filtfilt(b,a,double(NSP{i}.Data{SyncIndex(i)})');
    else
        NSP{i}.Data = filtfilt(b,a,double(NSP{i}.Data)');
    end
end

%Perform Cross Correlation
lagdiff = [];
SamplingRate = 30000;
[acor,lag] = xcorr(NSP{1}.Data(1:30000), NSP{2}.Data(1:30000));
[~, I] = max(abs(acor));
firstlagdiff = lag(I)

for i = 1:NumberOfFiles
    if Cell == 1;
        DataTimestamps(i) = NSP{i}.MetaTags.Timestamp(SyncIndex(i));
    else
        DataTimestamps(i) = NSP{i}.MetaTags.Timestamp;
    end
end

% If this method correctly predicts that actual offset of the
% files then we have a winner
if abs(firstlagdiff) < (abs(max(DataTimestamps) - min(DataTimestamps))+10)
    disp('Low frequency component found.')
    Method = 1;
else
    disp('No Low frequency component. Moving to next method.')
    Method = 0;
end

if ForceLFMethod == 1
    Method = 1;
end

if ForcePredictionMethod == 1
    Method = 0;
    disp('Prediction Method Forced. Moving to Prediction Method.');
end


if Method == 1
    fprintf(ReportFID,'Mode:');
    fprintf(ReportFID,'Low Frequency');
    fprintf(ReportFID,'\n');
    for i = 1:NumberOfFiles
        DataLength(i) = length(NSP{i}.Data);
    end
    EndPoint = min(DataLength);
    clear DataLength

    %Perform Cross Correlation
    lagdiff = [];
    SamplingRate = 30000;
    [acor,lag] = xcorr(NSP{1}.Data(EndPoint-30000:EndPoint), NSP{2}.Data(EndPoint-30000:EndPoint));
    [~, I] = max(abs(acor));
    plot(acor);
    lagdiff = lag(I)
    clear SamplingRate
    clear acor
    clear lag
    clear I
  

    %If Lag Amount is Negative, NSP1 Data needs to be resampled at a higher
    %rate. If Positive, NSP2 needs to be resampled at a higher rate.
    if (lagdiff < 0)
        DataToResample = 1;
    elseif (lagdiff > 0)
        DataToResample = 2;
    else
        disp('Your drift amount appears to be 0; these data may be aligned. Please check manually')
        DataToResample = 0;
        return
    end
    
    
    disp('Calculations complete. Opening full data file for drift correction. This may take a long while.')
    fprintf(ReportFID,'lagdiff:');
    fprintf(ReportFID,num2str(lagdiff));
    fprintf(ReportFID,'\n');
    for i = 1:NumberOfFiles
        NSP{i} = openNSx(Filename{i});
        if Cell == 1;     
            NSP{i}.Data = NSP{i}.Data{SyncIndex(i)};
            NSP{i}.MetaTags.Timestamp = NSP{i}.MetaTags.Timestamp(SyncIndex(i));
        else

        end
    end
    DriftAmount = abs(lagdiff - firstlagdiff);
    fprintf(ReportFID,'Drift Amount:');
    fprintf(ReportFID,num2str(DriftAmount));
    fprintf(ReportFID,'\n');
    OriginalLength = length(NSP{DataToResample}.Data);
    ResamplePeriod = round(EndPoint/DriftAmount);
    % Our frequency adjustments are too small for the Matlab resampling to
    % handle
%     OriginalSamplingRate = 30000;
%     NewSamplingRate = 30000 + 30000/ResamplePeriod;
%     [p,q] = rat(NewSamplingRate/OriginalSamplingRate);
%     DataForResampling = resample(double(DataForResampling),p,q);
    TotalPeriods = round(EndPoint/ResamplePeriod); %This is basically just equal to drift amount, but it looks nice
    
    
   if TestMethod == 1
        tic
            RepeatingArray = ones(1,length(NSP{DataToResample}.Data));
            RepeatingArray(1:ResamplePeriod:end) = 2;
            NSP{DataToResample}.Data = repelem(NSP{DataToResample}.Data,1,RepeatingArray);
        toc
   else
        for period = 1:TotalPeriods
            tic
            NSP{DataToResample}.Data = [NSP{DataToResample}.Data(:,1:ResamplePeriod*period) NSP{DataToResample}.Data(:,ResamplePeriod*period:end)];
            disp([num2str(period) 'of' num2str(TotalPeriods)]);
            t = toc;
            EstimatedTimeLeft = t*(TotalPeriods-period);
            disp(strcat('Estimated Time Left: ',num2str(EstimatedTimeLeft),' Seconds'))
        end
   end
    
    
    %Correct disparity between the files caused by resynchronization
    for i = 1:NumberOfFiles
        NSP{i}.Data = [zeros(size(NSP{i}.Data,1),NSP{i}.MetaTags.Timestamp) NSP{i}.Data];
    end
    
    MethodComment = 'Sync Method: Low Frequency Method';
    NSP{1}.MetaTags.Comment(end-length(MethodComment)+1:end) = MethodComment;
    NSP{2}.MetaTags.Comment(end-length(MethodComment)+1:end) = MethodComment;
    
    saveNSxSync(NSP{1});
    saveNSxSync(NSP{2});
    
    
    %Handle Additional File Types
    for i = 1:NumberOfFiles
        ValidFileTypes = {};
            for idx = 1:NumberOfFiles
                FileExtTypes = {'.ns1' '.ns2' '.ns3' '.ns4' '.ns5' '.ns6'};
                ValidTypesIndices = [];
                for Type = 1:length(FileExtTypes)
                    if ~strcmpi(FileExtTypes{Type},NSP{1}.MetaTags.FileExt)
                        ValidTypesIndices = [ValidTypesIndices Type];
                    end
                end
                ValidFileTypes{idx} = ValidTypesIndices;
            end
            
            
        for Type = ValidFileTypes{i}
            TempFilename = fullfile(NSP{i}.MetaTags.FilePath, [NSP{i}.MetaTags.Filename FileExtTypes{Type}]);
            if exist(TempFilename)
                disp(strcat('Found:',TempFilename));
                TempStructure = openNSx(TempFilename,'noread');
                EndPacket = TempStructure.MetaTags.DataPoints;
                TimestampScale = 30000/TempStructure.MetaTags.SamplingFreq;
                clear TempStructure
                NSx = openNSx(TempFilename);
                if DataToResample == i
                    if and(and(iscell(NSx.Data),Cell==1),length(NSx.Data)==2)
                        %disp('This function does not work on paused data currently');

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
                        TempResamplePeriod = round((ResamplePeriod/TimestampScale)*TimestampScale);
                        TempEndPoint = round(EndPoint/TimestampScale);
                        TotalPeriods = floor(TempEndPoint/TempResamplePeriod);
                        if TotalPeriods > 1
                            for period = 1:TotalPeriods
                                NSx.Data = [NSx.Data(:,1:TempResamplePeriod*period) NSx.Data(:,TempResamplePeriod*period:end)];
                                disp([num2str(period) 'of' num2str(TotalPeriods)]);
                            end
                        end
                    end
                end
                
                
                saveNSxSync(NSx)
            end
        end
        
        
        %NEV Attempt
        if (exist(fullfile(NSP{i}.MetaTags.FilePath, [NSP{i}.MetaTags.Filename '.nev'])))
            TempNEV = openNEV(fullfile(NSP{i}.MetaTags.FilePath, [NSP{i}.MetaTags.Filename '.nev']),'nosave','nomat');
            disp(strcat('Found:',fullfile(NSP{i}.MetaTags.FilePath, [NSP{i}.MetaTags.Filename '.nev'])));
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

            
            saveNEVSync(TempNEV)
        end
    end
    
    
    fclose(ReportFID);
    return
else
    %disp('No low frequency component moving on to prediction method.')
end 


%% Complete Prediction Method


%Top PFC
%Bottom FEF (Analog Inputs)

fprintf(ReportFID,'Mode:');
fprintf(ReportFID,'Prediction');
fprintf(ReportFID,'\n');

for i = 1:NumberOfFiles
    if Cell == 1
        DataTimestamps(i) = NSP{i}.MetaTags.Timestamp(SyncIndex(i));
    else
        DataTimestamps(i) = NSP{i}.MetaTags.Timestamp;
    end
end

%Force method to 1 because this is the last option
Method = 1;



if Method == 1 
    for i = 1:NumberOfFiles
        DataLength(i) = length(NSP{i}.Data);
    end
    EndPoint = min(DataLength);
    clear DataLength
    
    
    AnalogInputCounts = [];
    for i = 1:NumberOfFiles
        AnalogInputCounts(i) = length(find([NSP{i}.ElectrodesInfo.ElectrodeID]>128));
    end
    
    %Check for event codes in digital inputs of NEV file. This will replace
    %the section below
    for i = 1:NumberOfFiles
        if (exist(fullfile(NSP{i}.MetaTags.FilePath, [NSP{i}.MetaTags.Filename '.nev'])))

        end
    end
    
    DataToResample = [];
    for i = 1:NumberOfFiles
        if (exist(fullfile(NSP{i}.MetaTags.FilePath, [NSP{i}.MetaTags.Filename '.nev'])))
            TempNEV = openNEV(fullfile(NSP{i}.MetaTags.FilePath, [NSP{i}.MetaTags.Filename '.nev']),'nosave','nomat');
            disp(strcat('Found:',fullfile(NSP{i}.MetaTags.FilePath, [NSP{i}.MetaTags.Filename '.nev'])));
            if ~isempty(NEV.Data.SerialDigitalIO.TimeStamp) 
                if isempty(DataToResample)
                    DataToResample = i;
                else
                    disp('Multiple files contain digital events. Cannot distinguish NSP 1.5');
                    return
                end
            end
        end 
    end
    
    clear TempNEV

    for i = 1:NumberOfFiles
        disp('Calculations complete. Opening full data file for drift correction. This may take a long while.')
        NSP{i} = openNSx(Filename{i});
        if Cell == 1
            NSP{i}.Data = NSP{i}.Data{SyncIndex(i)};
            NSP{i}.MetaTags.Timestamp = NSP{i}.MetaTags.Timestamp(SyncIndex(i));
        end
    end
    
    OriginalLength = length(NSP{DataToResample}.Data);
    DriftAmount = 100*(OriginalLength/(30000*60*5)); %We see around 100 samples loss every 5 minutes on most of these NSPs.
    fprintf(ReportFID,'Drift Amount:');
    fprintf(ReportFID,num2str(DriftAmount));
    fprintf(ReportFID,'\n');
    ResamplePeriod = round(EndPoint/DriftAmount);
    
    TotalPeriods = EndPoint/ResamplePeriod; %This is basically just equal to drift amount, but it looks nice

    if TestMethod == 1
        tic
            RepeatingArray = ones(1,length(NSP{DataToResample}.Data));
            RepeatingArray(1:ResamplePeriod:end) = 2;
            NSP{DataToResample}.Data = repelem(NSP{DataToResample}.Data,1,RepeatingArray);
        toc
    else
        for period = 1:TotalPeriods
            tic
            NSP{DataToResample}.Data = [NSP{DataToResample}.Data(:,1:ResamplePeriod*period) NSP{DataToResample}.Data(:,ResamplePeriod*period:end)];
            disp([num2str(period) 'of' num2str(TotalPeriods)]);
            t = toc;
            EstimatedTimeLeft = t*(TotalPeriods-period);
            disp(strcat('Estimated Time Left: ',num2str(EstimatedTimeLeft),' Seconds'))
        end
    end


    
    %Correct disparity between the files caused by resynchronization
    for i = 1:NumberOfFiles
        NSP{i}.Data = [zeros(size(NSP{i}.Data,1),NSP{i}.MetaTags.Timestamp) NSP{i}.Data];
    end
    
    MethodComment = 'Sync Method: Prediction Method';
    NSP{1}.MetaTags.Comment(end-length(MethodComment)+1:end) = MethodComment;
    NSP{2}.MetaTags.Comment(end-length(MethodComment)+1:end) = MethodComment;
    
    
    saveNSxSync(NSP{1});
    saveNSxSync(NSP{2});
    
    %Handle Additional File Types
    for i = 1:NumberOfFiles
        ValidFileTypes = {};
            for idx = 1:NumberOfFiles
                FileExtTypes = {'.ns1' '.ns2' '.ns3' '.ns4' '.ns5' '.ns6'};
                ValidTypesIndices = [];
                for Type = 1:length(FileExtTypes)
                    if ~strcmpi(FileExtTypes{Type},NSP{1}.MetaTags.FileExt)
                        ValidTypesIndices = [ValidTypesIndices Type];
                    end
                end
                ValidFileTypes{idx} = ValidTypesIndices;
            end
            
            
        for Type = ValidFileTypes{i}
            TempFilename = fullfile(NSP{i}.MetaTags.FilePath, [NSP{i}.MetaTags.Filename FileExtTypes{Type}]);
            if exist(TempFilename)
                disp(strcat('Found:',TempFilename));
                TempStructure = openNSx(TempFilename,'noread');
                EndPacket = TempStructure.MetaTags.DataPoints;
                TimestampScale = 30000/TempStructure.MetaTags.SamplingFreq;
                clear TempStructure
                NSx = openNSx(TempFilename);
                if DataToResample == i
                    if and(and(iscell(NSx.Data),Cell==1),length(NSx.Data)==2)
                        %disp('This function does not work on paused data currently');

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
                        TempResamplePeriod = round((ResamplePeriod/TimestampScale)*TimestampScale);
                        TempEndPoint = round(EndPoint/TimestampScale);
                        TotalPeriods = floor(TempEndPoint/TempResamplePeriod);
                        if TotalPeriods > 1
                            for period = 1:TotalPeriods
                                NSx.Data = [NSx.Data(:,1:TempResamplePeriod*period) NSx.Data(:,TempResamplePeriod*period:end)];
                                disp([num2str(period) 'of' num2str(TotalPeriods)]);
                            end
                        end
                    end
                end
                
                
                saveNSxSync(NSx)
            end
        end
        
        
        %NEV Attempt
        if (exist(fullfile(NSP{i}.MetaTags.FilePath, [NSP{i}.MetaTags.Filename '.nev'])))
            TempNEV = openNEV(fullfile(NSP{i}.MetaTags.FilePath, [NSP{i}.MetaTags.Filename '.nev']),'nosave','nomat');
            disp(strcat('Found:',fullfile(NSP{i}.MetaTags.FilePath, [NSP{i}.MetaTags.Filename '.nev'])));
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

            
            saveNEVSync(TempNEV)
        end
    end
    
    
    fclose(ReportFID);
    return
else
    disp('No other options. No Data was synced.')
end 

%% Wrap Up
