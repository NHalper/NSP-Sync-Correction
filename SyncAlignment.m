function SyncAlignment(NumberOfFiles)

%%% Alignment Script for Multiple NSPs
%
% Purpose: This script produces files that are modified so that alignment
% of data from multiple NSPs is intuitive to understand. It does this by
% ensuring that the timestamps are on the same scale between all files.
%
% Author: Nick Halper for Blackrock Microsystems
% Contact: nhalper@blackrockmicro.com
% www.blackrockmicro.com
%
% Version 1.0
% V1.00 = Initial Release
% V1.01 = Corrected issue that would allow the function to read non-30k files.
% V1.02 = Corrected an issue that would have allowed NS6 files to incorrectly erase their non-zero timestamp starts if an NS5 file was used as the base file


NSPInfo = {};
uiwait(msgbox('Please select the files to align. One at a time.', 'Info','modal'));
for i = 1 :NumberOfFiles
    NSPInfo{i} = openNSx('noread');
    if ~isstruct(NSPInfo{i})
        disp('No File Selected. Terminating Script.')
        return
    end
    if length(NSPInfo{i}.MetaTags.Timestamp) > 1
        disp('This function does not deal with paused data or missing segments, please contact support')
        return
    end
    if NSPInfo{i}.MetaTags.SamplingFreq ~= 30000
        disp('This function requires inputs of 30k sampled data.')
        return
    end
end

NSPFilenames = {};
GenericNSPFilenames = {};
for i = 1:NumberOfFiles
    NSPFilenames{i} = fullfile(NSPInfo{i}.MetaTags.FilePath,[NSPInfo{i}.MetaTags.Filename NSPInfo{i}.MetaTags.FileExt]);
    GenericNSPFilenames{i} = NSPFilenames{i}(1:end-4);
end

%Get sync Pulse Data
for i = 1:NumberOfFiles
    NSPCodeData{i} = openNSx(NSPFilenames{i},['c:' num2str(NSPInfo{i}.MetaTags.ChannelCount)]);
    NSPCodeData{i}.Data = [repmat(0,size(NSPCodeData{i}.Data,1),NSPInfo{i}.MetaTags.Timestamp) NSPCodeData{i}.Data];
    NSPCodes{i} = syncPatternDetectNSx(NSPCodeData{i}.Data(size(NSPCodeData{i}.Data,1),1:30*NSPInfo{i}.MetaTags.SamplingFreq));
end
clear NSPCodeData


%Find timestamp of earliest common code
Success = 0;
FinalIntersection = [];

TempIntersection = NSPCodes{1}{1};
for i = 2:NumberOfFiles
    TempIntersection = intersect(TempIntersection,NSPCodes{i}{1});
end
FinalIntersection = TempIntersection;
clear TempIntersection;
FirstCommonCode = FinalIntersection(1);
for i = 1:NumberOfFiles
    CodeTimestamps(i) = NSPCodes{i}{2}(find(NSPCodes{i}{1}==FirstCommonCode,1));
end
    
%Create array of additional file types to check
FileExtTypes = {'.ns1' '.ns2' '.ns3' '.ns4' '.ns5' '.ns6'};
ValidTypesIndices = [];
for Type = 1:length(FileExtTypes)
    if ~strcmpi(FileExtTypes{Type},NSPInfo{1}.MetaTags.FileExt)
        ValidTypesIndices = [ValidTypesIndices Type];
    end
end


%Perform alignment and save new files; could remove the above File
%Extension Array and merge the leading lines with the For loop in this
%array
for i = 1:NumberOfFiles
    NSPData = openNSx(NSPFilenames{i});
    TimestampScale = 30000/NSPData.MetaTags.SamplingFreq;
    %NSPData.Data = [repmat(0,size(NSPData.Data,1),(max(CodeTimestamps)-CodeTimestamps(i))+NSPData.MetaTags.Timestamp) NSPData.Data];
    NSPData.Data = [repmat(0,size(NSPData.Data,1),round((max(CodeTimestamps)-CodeTimestamps(i))/TimestampScale + NSPData.MetaTags.Timestamp)) NSPData.Data];
    NSPData.MetaTags.Timestamp = 0;
    saveNSxSync(NSPData);
    for Type = ValidTypesIndices
        if exist([GenericNSPFilenames{i} FileExtTypes{Type}])
            NSPData = openNSx([GenericNSPFilenames{i} FileExtTypes{Type}]);
            TimestampScale = 30000/NSPData.MetaTags.SamplingFreq;
            NSPData.Data = [repmat(0,size(NSPData.Data,1),round((max(CodeTimestamps)-CodeTimestamps(i))/TimestampScale)) NSPData.Data];
            NSPData.MetaTags.Timestamp = 0;
            saveNSxSync(NSPData);
        end
    end
end

    

%Repeat process for NEV Files
for i = 1:NumberOfFiles
    if (exist([GenericNSPFilenames{i} '.nev']))
        TempNEV = openNEV([GenericNSPFilenames{i} '.nev'],'nosave','nomat');

        TempNEV.Data.SerialDigitalIO.TimeStamp = TempNEV.Data.SerialDigitalIO.TimeStamp + (max(CodeTimestamps)-CodeTimestamps(i));
        TempNEV.Data.Spikes.TimeStamp = TempNEV.Data.Spikes.TimeStamp + (max(CodeTimestamps)-CodeTimestamps(i));
        TempNEV.Data.Comments.TimeStamp = TempNEV.Data.Comments.TimeStamp + (max(CodeTimestamps)-CodeTimestamps(i));
        TempNEV.Data.VideoSync.TimeStamp = TempNEV.Data.VideoSync.TimeStamp + (max(CodeTimestamps)-CodeTimestamps(i));
        if ~isempty(TempNEV.Data.Tracking)
            uiwait(msgbox('This version does not operate on tracking data timestamps. Please contact Blackrock Support.','ERROR','modal'));
        end
        TempNEV.Data.TrackingEvents.TimeStamp = TempNEV.Data.TrackingEvents.TimeStamp + (max(CodeTimestamps)-CodeTimestamps(i));
        TempNEV.Data.PatientTrigger.TimeStamp = TempNEV.Data.PatientTrigger.TimeStamp + (max(CodeTimestamps)-CodeTimestamps(i));
        TempNEV.Data.Reconfig.TimeStamp = TempNEV.Data.Reconfig.TimeStamp + (max(CodeTimestamps)-CodeTimestamps(i));

        saveNEVSync(TempNEV)
    end
end

