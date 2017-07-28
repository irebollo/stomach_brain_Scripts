function [ markersDownsampled ] = prepro_egg_downsampleMarkers(data_raw,data_downsampled,markers_raw)

%{
Get closest sample number of new marker based on the timestamps of the not downsampled data

INPUTS:

data_raw - data not downsampled
data_downsampled - downsampled data
markers_raw - the original event structure

OUTPUTS:
markersDownsampled: the event structure with the synchronized timestamps

Commented IR 27/06/2017

%}

fprintf('\n###############\nSynchro dsp events...\n\n')


    %search for timestamp of each sample for each event
    event_start = markers_raw.trl(:,1);
    event_end = markers_raw.trl(:,2);
    markers_raw.trl(:,3) = data_raw.time{1,1}(1,event_start);
    markers_raw.trl(:,4) = data_raw.time{1,1}(1,event_end);
    
    markersDownsampled = markers_raw.trl;
    
    % for each event in the eventstructure seearch it's corresponding
    % timestamp in the downsampled data
    for iEvent=1:length(markers_raw.trl)
        index_newBegin  = nearest(data_downsampled.time{1},double(markers_raw.trl(iEvent,3)));
        index_newEnd  = nearest(data_downsampled.time{1},double(markers_raw.trl(iEvent,4)));
        markersDownsampled(iEvent,1) = data_downsampled.time{1}(index_newBegin);
        markersDownsampled(iEvent,2) = data_downsampled.time{1}(index_newEnd);
        markersDownsampled(iEvent,3) = index_newBegin;
        markersDownsampled(iEvent,4) = index_newEnd;


    end
fprintf('\n###############\nDONE\n\n')
end
