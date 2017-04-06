function splitFileResets

% splitFileResets
% 
% Opens and splits both an NEV and NSx file in pieces timewise at clock
% restarts.
%
% Use splitFileResets

% All input arguments are optional. Input arguments can be in any order.
%
%
%   Example 1: 
%   splitFileResets;
%
%   In the example above, the files will be split at every clock restart in
%   the recording.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Nick Halper
%   support@blackrockmicro.com
%   Blackrock Microsystems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Version History
%
% 1.0.0.0:
%   - Initial release.
%
% 1.1.0.0:
%   - Improved RAM usage and corrected a few bugs with segment assigments
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

% Validating input parameter
% if ~exist('splitCount', 'var')
%     splitCount = 2;
% end

% Getting the file name
[contfname, contpath] = getFile('*.NS6', 'Choose an NSx file...');
if contfname == 0
    disp('No file was selected.');
    if nargout
        clear variables;
    end
    return;
end
contfext = contfname(end-3:end);

% Getting the file name
[fname, path] = getFile('*.nev', 'Choose an NEV file...');
if fname == 0
    disp('No file was selected.');
    if nargout
        clear variables;
    end
    return;
end
fext = fname(end-3:end);


% Loading the file
%% Reading Basic Header from file.

%Reading NSx header
contFID                     = fopen([contpath contfname], 'r', 'ieee-le');
contBasicHeader             = fread(contFID, 314, '*uint8');
contpositionEOE             = typecast(contBasicHeader(11:14), 'uint32');
contchannelCount            = double(typecast(contBasicHeader(311:314), 'uint32'));
fseek(contFID,0,'bof')
contTotalHeaders            = fread(contFID,contpositionEOE,'char');
fseek(contFID,0,'eof')
contEndOfFile               = ftell(contFID);
fseek(contFID,contpositionEOE,'bof')

%Reading NEV header
FID                       = fopen([path fname], 'r', 'ieee-le');
BasicHeader               = fread(FID, 336, '*uint8');
Trackers.fExtendedHeader  = double(typecast(BasicHeader(13:16), 'uint32'));
Trackers.countPacketBytes = double(typecast(BasicHeader(17:20), 'uint32'));

%% Doing Trackers

%NSx Trackers
contSegmentCount = 0;
contDataPointsInSegment = [];
contTimestamps = [];
RestartIndex = [1];
while ftell(contFID)<contEndOfFile
    try
    contSegmentCount = contSegmentCount+1;
    %Read data header
    contDataHeader = fread(contFID,9,'char');
    fseek(contFID,-9,'cof');
    fseek(contFID,1,'cof');
    contTimestamps(contSegmentCount) = fread(contFID,1,'uint32');
    contNumberOfDataPoints = fread(contFID,1,'uint32');
    contDataPointsInSegment(contSegmentCount) = contNumberOfDataPoints;
    %Read data
    fseek(contFID,contNumberOfDataPoints*contchannelCount* 2,'cof');
    ftell(contFID)
    if ftell(contFID) == contEndOfFile
        disp('done')
    end
    if contSegmentCount > 1
        if contTimestamps(contSegmentCount-1)+contDataPointsInSegment(contSegmentCount-1)>contTimestamps(contSegmentCount)
            RestartIndex = [RestartIndex contSegmentCount];
        end
    end
    catch
        disp('Error in counting data segments');
    end
end
fseek(contFID,contpositionEOE,'bof')
fprintf('Found %d Segments in Continuous Data File. ', length(RestartIndex));


%NEV Trackers
fseek(FID, 0, 'eof');
Trackers.fData = ftell(FID);
Trackers.countDataPacket = (Trackers.fData - Trackers.fExtendedHeader)/Trackers.countPacketBytes;


fseek(FID, Trackers.fExtendedHeader, 'bof');
tRawData  = fread(FID, [10 Trackers.countDataPacket], '10*uint8=>uint8', Trackers.countPacketBytes - 10);
Timestamp = tRawData(1:4,:);
Timestamp = typecast(Timestamp(:), 'uint32').';

splitPacketStarts = [find(diff(double(Timestamp))<0) length(Timestamp)];
%splitPacketBytes = Trackers.countPacketBytes * splitPacketStarts;
splitPacketBytes = Trackers.countPacketBytes * [splitPacketStarts(1) diff(splitPacketStarts)];

fprintf('Found %d Segments in NEV Data File. ', length(splitPacketStarts));

%% Comparing Files and Matching Segments
if length(splitPacketStarts) == contSegmentCount
    disp('Segment counts match. Continuing to file writing...')
else
    disp('Segment counts do not match. Attempting to match segments...')
    
    contTotalTimestamps = [];
    for idx = 1:contSegmentCount
        if ~isempty(find(RestartIndex==idx))
           
            contCurrentRestartIndex = find(RestartIndex==idx);
            contStartingTimestamp = contTimestamps(idx);
        end
        if (~isempty(find(RestartIndex==idx+1))) || (idx == contSegmentCount)
            contTotalTimestamps(contCurrentRestartIndex) = contTimestamps(idx)+contDataPointsInSegment(idx);
        end
    end
    disp('Continuous Data Segment Lengths:')
    disp(contTotalTimestamps)
    nevTotalTimestamps = Timestamp([splitPacketStarts(2:end) length(Timestamp)]);
    disp('NEV Data Segment Lengths:')
    disp(nevTotalTimestamps)
    
    nevPairedContinuousSegment = [];
    tempcontSegmentCounter = 1;
    for idx = 1:length(nevTotalTimestamps)
        DoNotRepeat = 0;
        for contSegment = tempcontSegmentCounter:length(contTotalTimestamps)
            if nevTotalTimestamps(idx) < contTotalTimestamps(contSegment) && DoNotRepeat == 0
                nevPairedContinuousSegment(idx) = contSegment;
                tempcontSegmentCounter = contSegment+1;
                DoNotRepeat = 1;
            end
        end
    end
    disp('Matches assigned.');
    %return
end



%% Writing Files

% Writing NSx File
for idx = 1:contSegmentCount
    if ~isempty(find(RestartIndex==idx))
        contFIDw                    = fopen([contpath contfname(1:end-4) '-s' sprintf('%03d',find(RestartIndex==idx)) contfname(end-3:end)], 'w+', 'ieee-le');
    end
    fprintf('\nReading NSx segment %d... ', idx);
    contdataHeader              = fread(contFID,9,'char');
    fseek(contFID,-9,'cof');
    fseek(contFID,5,'cof');
    contNumberOfDataPoints      = fread(contFID,1,'uint32');
    contdataSegment             = fread(contFID,contNumberOfDataPoints*contchannelCount* 2,'*uint8');
    fprintf('Writing NSx segment %d... ', idx);
    if ~isempty(find(RestartIndex==idx))
        fwrite(contFIDw, contTotalHeaders, 'char');
    end
    if not(contDataPointsInSegment(idx) == 0)
        fwrite(contFIDw, contdataHeader, 'char');
        fwrite(contFIDw, contdataSegment, '*uint8');
    end
    clear contdataSegment;
    if (~isempty(find(RestartIndex==idx+1))) || (idx == contSegmentCount)
        fclose(contFIDw);
    end
end

% Writing NEV File
% Reading headers and seeking to beginning of data
fseek(FID, 0, 'bof');
fileHeader = fread(FID, Trackers.fExtendedHeader, '*uint8');

for idx = 1:length(splitPacketStarts)
    % Opening a file for saving
    FIDw = fopen([path fname(1:end-4) '-s' sprintf('%03d', nevPairedContinuousSegment(idx)) fname(end-3:end)], 'w+', 'ieee-le');
    fprintf('\nReading NEV segment %d... ', idx);
    % Reading the segment
    dataSegment = fread(FID, splitPacketBytes(idx), 'char');
    fprintf('Writing NEV segment %d... ', idx);
    % Writing the segmented data into file
    fwrite(FIDw, fileHeader, 'char');
    fwrite(FIDw, dataSegment, 'char');
    % Clearing variables and closing file
    clear dataSegment;
    fclose(FIDw);
end
