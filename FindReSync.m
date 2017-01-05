function ResyncIndices = FindReSync(NSxStructureMetaTags)
ResyncIndices = [];
if length(NSxStructureMetaTags.Timestamp)>1
    for idx = 1:length(NSxStructureMetaTags.Timestamp)-1
        if NSxStructureMetaTags.Timestamp(idx)+NSxStructureMetaTags.DataPoints(idx)*NSxStructureMetaTags.SamplingFreq/30000>NSxStructureMetaTags.Timestamp(idx+1)
            ResyncIndices = [ResyncIndices idx+1];
        end
    end
end

if isempty(ResyncIndices)
    disp('No resync events detected. This file may have pauses or packet loss.')
    ResyncIndices = 1;
end

if length(ResyncIndices)>1
    disp('Multiple resync events found. This could indicate a problem with the data file, please inspect manually.')
    ResyncIndices = 0;
end




