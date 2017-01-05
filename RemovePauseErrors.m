function NSxStructure = RemovePauseErrors(NSxStructure)

SyncIndex = FindReSync(NSxStructure.MetaTags);
if SyncIndex == 1
    %File has pauses or packet loss
    %Proceed as normal
elseif SyncIndex == 0
    %File has multiple resyncs. This could indicate an issue with file
    return
else
    %Use Sync Index to Clean the Data
    NSxStructure.MetaTags.Timestamp = NSxStructure.MetaTags.Timestamp(SyncIndex:end);
    NSxStructure.MetaTags.DataPoints = NSxStructure.MetaTags.Timestamp(SyncIndex:end);
    NSxStructure.Data = NSxStructure.Data(SyncIndex:end);
end

%Find and remove pause data
BadSectionsIndices = find(NSxStructure.MetaTags.DataPoints<2);

NSxStructure.MetaTags.Timestamp(BadSectionsIndices) = [];
NSxStructure.MetaTags.DataPoints(BadSectionsIndices) = [];
NSxStructure.Data(BadSectionsIndices) = [];

saveNSxSync(NSxStructure);

return

end