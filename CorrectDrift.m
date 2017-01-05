function CorrectDrift(NumberOfFiles)

%%% Drift Correction Script for Jake and Roman
%Version 1.02
%Change Log
%V1.00 = Initial release
%V1.01 = Nearly complete rewrite; Rework to allow multiple NSPs
%V1.02 = Improved NEV functionality
%V1.03 = Better ability to handle 256 channel synchronized systems.
%V1.04 = Drastically improved speed. Up to 50x faster.
NSPSyncSignalStructures = {};
uiwait(msgbox('Please select the files to align. One at a time.','Info','modal'));
for i = 1:NumberOfFiles
    NSPSyncSignalStructures{i} = openNSx;
    if NSPSyncSignalStructures{i}.MetaTags.SamplingFreq ~= 30000
        disp('Sync pulse is currently required to be sampled at 30kHz')
        return
    end
end

for i = 1:NumberOfFiles
    if iscell(NSPSyncSignalStructures{i}.Data)
        disp('Cells in data');
        SyncIndex = FindReSync(NSPSyncSignalStructures{i}.MetaTags);
        if SyncIndex == 0
            return
        end
        if SyncIndex < length(NSPSyncSignalStructures{i}.Data)-1
            disp('File may have packet loss. This version of the script does not currently work with packet loss.')
            return
        end
        if SyncIndex == length(NSPSyncSignalStructures{i}.Data)
            NSPSyncSignalStructures{i}.Data = [zeros(size(NSPSyncSignalStructures{i}.Data{SyncIndex},1),NSPSyncSignalStructures{i}.MetaTags.Timestamp(SyncIndex)) NSPSyncSignalStructures{i}.Data{SyncIndex}];
            NSPSyncSignalStructures{i}.MetaTags.Timestamp = 0;
            NSPSyncSignalStructures{i}.MetaTags.DataPoints = length(NSPSyncSignalStructures{i}.Data);
        end
    else
        NSPSyncSignalStructures{i}.Data = [zeros(size(NSPSyncSignalStructures{i}.Data,1),NSPSyncSignalStructures{i}.MetaTags.Timestamp) NSPSyncSignalStructures{i}.Data];
        NSPSyncSignalStructures{i}.MetaTags.Timestamp = 0;
        NSPSyncSignalStructures{i}.MetaTags.DataPoints = length(NSPSyncSignalStructures{i}.Data);
    end
end



NSPFilenames = {};
for i = 1:NumberOfFiles
    NSPFilenames{i} = fullfile(NSPSyncSignalStructures{i}.MetaTags.FilePath,[NSPSyncSignalStructures{i}.MetaTags.Filename NSPSyncSignalStructures{i}.MetaTags.FileExt]);
end

GenericNSPFilenames = {};
for i = 1:NumberOfFiles
    GenericNSPFilenames{i} = NSPFilenames{i}(1:end-4);
end


%Find the latest point at which we can calculate drift
FileLengths = [];
for i = 1:NumberOfFiles
    FileLengths = [FileLengths length(NSPSyncSignalStructures{i}.Data)];
end
EndPoint = min(FileLengths);

%Find codes in sync pulse
for i = 1:NumberOfFiles
    CodeResults{i} = syncPatternDetectNSx(NSPSyncSignalStructures{i}.Data(end,:));
end

%Find First Common Code
IntersectResult = intersect(CodeResults{1}{1}(1:3),CodeResults{2}{1}(1:3));
if NumberOfFiles>2  
    for i = 3:NumberOfFiles
        IntersectResult = intersect(IntersectResult,CodeResults{i}{1}(1:3));
    end
end
FirstCommonCode = IntersectResult(1);
FirstCodeTimestamps = [];
for i = 1:NumberOfFiles
    CodeIndex = find(CodeResults{i}{1}(1:3)==FirstCommonCode);
    FirstCodeTimestamps(i) = CodeResults{i}{2}(CodeIndex);
end
    

%Find Last Common Code
CodeLengths = [];
for i = 1:NumberOfFiles
    CodeLengths = [CodeLengths length(CodeResults{i}{1})];
end
FinalCodeIndex = min(CodeLengths);

IntersectResult = intersect(CodeResults{1}{1}(FinalCodeIndex-3:FinalCodeIndex),CodeResults{2}{1}(FinalCodeIndex-3:FinalCodeIndex));
if NumberOfFiles>2  
    for i = 3:NumberOfFiles
        IntersectResult = intersect(IntersectResult,CodeResults{i}{1}(FinalCodeIndex-3:FinalCodeIndex));
    end
end
LastCommonCode = max(IntersectResult);
LastCodeTimestamps = [];
for i = 1:NumberOfFiles
    CodeIndex = find(CodeResults{i}{1}(FinalCodeIndex-3:FinalCodeIndex)==LastCommonCode);
    LastCodeTimestamps(i) = CodeResults{i}{2}(CodeIndex+FinalCodeIndex-4);
end
    
%Check that data start out aligned. 
if(max(FirstCodeTimestamps)-min(FirstCodeTimestamps))>6
    %First code may be late enough that some clock drift has occurred. Set
    %to 6 to account for this, but it should be rare.
    disp('Original data does not seem to be aligned. Please use alignment script beforehand.')
    return
end

%Find Max Length Data. Resample other data to meet Max Length Data
LatestCodeTimestamp = max(LastCodeTimestamps);
DriftAmount = [];
for i = 1:NumberOfFiles
    DriftAmount(i) = LatestCodeTimestamp-LastCodeTimestamps(i);
end

if sum(DriftAmount) == 0
    disp('Your drift amount appears to be zero. These data files may not have drift. Inspect manually.')
end


%Correct for drift and propgate to other files.
for i = 1:NumberOfFiles
    if DriftAmount(i) ~= 0
        OriginalLength = length(NSPSyncSignalStructures{i}.Data);
        ResamplePeriod = round(EndPoint/DriftAmount(i));
        OriginalSamplingRate = NSPSyncSignalStructures{i}.MetaTags.SamplingFreq; %Should be 30k for now.
        NewSamplingRate = OriginalSamplingRate + OriginalSamplingRate/ResamplePeriod;
        [p,q] = rat(NewSamplingRate/OriginalSamplingRate);
%         if p*q <  2^31   
%             NSPSyncSignalStructures{i}.Data = resample(double(NSPSyncSignalStructures{i}.Data),p,q);
%         else

            RepeatingArray = ones(1,length(NSPSyncSignalStructures{i}.Data));
            RepeatingArray(1:ResamplePeriod:end) = 2;
            NSPSyncSignalStructures{i}.Data = repelem(NSPSyncSignalStructures{i}.Data,1,RepeatingArray);

            %Old Method = Slow
%             TotalPeriods = EndPoint/ResamplePeriod;
%             for period = 1:TotalPeriods
%                 NSPSyncSignalStructures{i}.Data = [NSPSyncSignalStructures{i}.Data(:,1:ResamplePeriod*period) NSPSyncSignalStructures{i}.Data(:,ResamplePeriod*period:end)];
%                 disp([num2str(period) 'of' num2str(TotalPeriods)]);
%             end
%         end
        if abs((length(NSPSyncSignalStructures{i}.Data) - OriginalLength)-DriftAmount) > 2
            disp('There may be an issue with the resampling of the data, please manually check it.')
           
        end
        
        ValidFileTypes = {};
        for idx = 1:NumberOfFiles
            FileExtTypes = {'.ns1' '.ns2' '.ns3' '.ns4' '.ns5' '.ns6'};
            ValidTypesIndices = [];
            for Type = 1:length(FileExtTypes)
                if ~strcmpi(FileExtTypes{Type},NSPSyncSignalStructures{idx}.MetaTags.FileExt)
                    ValidTypesIndices = [ValidTypesIndices Type];
                end
            end
            ValidFileTypes{idx} = ValidTypesIndices;
        end
        
        for Type = ValidFileTypes{i}
            if (exist([GenericNSPFilenames{i} FileExtTypes{Type}]))
                TempStructure = openNSx([GenericNSPFilenames{i} FileExtTypes{Type}],'noread');
                EndPacket = TempStructure.MetaTags.DataPoints;
                TimestampScale = 30000/TempStructure.MetaTags.SamplingFreq;
                clear TempStructure
                NSx = openNSx([GenericNSPFilenames{i} FileExtTypes{Type}]);
                if iscell(NSx.Data)
                    %disp('This function does not work on paused data currently');
                    
                    NSx.Data = [zeros(size(NSx.Data{end},1),round(NSx.MetaTags.Timestamp(end)/30000)) NSx.Data{end}];
                    NSx.MetaTags.Timestamp = 0;
                    TempResamplePeriod = round((ResamplePeriod/TimestampScale)*TimestampScale);
                    TempEndPoint = round(EndPoint/TimestampScale);
                    TotalPeriods = floor(TempEndPoint/TempResamplePeriod);
                    if TotalPeriods > 1
                        
                        RepeatingArray = ones(1,length(NSx.Data));
                        RepeatingArray(1:TempResamplePeriod:end) = 2;
                        NSx.Data = repelem(NSx.Data,1,RepeatingArray); 
%                         Old Method is Slow
%                         for period = 1:TotalPeriods
%                             NSx.Data = [NSx.Data(:,1:TempResamplePeriod*period) NSx.Data(:,TempResamplePeriod*period:end)];
%                             disp([num2str(period) 'of' num2str(TotalPeriods)]);
%                         end
                        
                    end
                    
                else
                    NSx.Data = [zeros(size(NSx.Data,1),NSx.MetaTags.Timestamp) NSx.Data];
                    TempResamplePeriod = round((ResamplePeriod/TimestampScale)*TimestampScale);
                    TempEndPoint = round(EndPoint/TimestampScale);
                    TotalPeriods = floor(TempEndPoint/TempResamplePeriod);
                    if TotalPeriods > 1
                        
                        RepeatingArray = ones(1,length(NSx.Data));
                        RepeatingArray(1:TempResamplePeriod:end) = 2;
                        NSx.Data = repelem(NSx.Data,1,RepeatingArray);  
%                         
%                         for period = 1:TotalPeriods
%                             NSx.Data = [NSx.Data(:,1:TempResamplePeriod*period) NSx.Data(:,TempResamplePeriod*period:end)];
%                             disp([num2str(period) 'of' num2str(TotalPeriods)]);
%                         end
                        
                    end
                end
                
                
                saveNSxSync(NSx)
            end
        end
        
        
        if (exist([GenericNSPFilenames{i} '.nev']))
            TempNEV = openNEV([GenericNSPFilenames{i} '.nev'],'nosave','nomat');
            
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
        
    else %Drift value is 0
        ValidFileTypes = {};
        for idx = 1:NumberOfFiles
            FileExtTypes = {'.ns1' '.ns2' '.ns3' '.ns4' '.ns5' '.ns6'};
            ValidTypesIndices = [];
            for Type = 1:length(FileExtTypes)
                if ~strcmpi(FileExtTypes{Type},NSPSyncSignalStructures{idx}.MetaTags.FileExt)
                    ValidTypesIndices = [ValidTypesIndices Type];
                end
            end
            ValidFileTypes{idx} = ValidTypesIndices;
        end
        
        %Save the file without modifying data, for consistency.
        for Type = ValidFileTypes{i}
            if (exist([GenericNSPFilenames{i} FileExtTypes{Type}]))
                NSx = openNSx([GenericNSPFilenames{i} FileExtTypes{Type}]);
                if iscell(NSx.Data)
                    NSx.Data = [zeros(size(NSx.Data{end},1),NSx.MetaTags.Timestamp(end)) NSx.Data{end}];
                    NSx.MetaTags.Timestamp = 0;
                    NSx.MetaTags.DataPoints = length(NSx.Data);
                end
                saveNSxSync(NSx)
            end
        end 
        if (exist([GenericNSPFilenames{i} '.nev']))
            TempNEV = openNEV([GenericNSPFilenames{i} '.nev'],'nosave','nomat');
            saveNEVSync(TempNEV)
        end
    end
    saveNSxSync(NSPSyncSignalStructures{i});
end
