function [ events_eyelink ] = GUTSEEMEG_prepare_eyelink_events(raw_data_path)
%GUTSEE_PREPARE_EYELINK_EVENTS Prepares eyelink event strucure
%   The original event structre obtained by calling my_read_eyelink_asc
%   does not contain timestamps for blinks and saccades, and does not note
%   triggers in the appropriate format (names and samples are noted in long
%   strings and thus difficult to access). This is corrected for in this
%   function. The resulting event structure is used for further analysis.
%
% Author: Nicolai Wolpert
% Email: Nicolai.Wolpert@ens.fr
% Version: 15.01.2018
% 
% INPUT:
% raw_data_path
% nBlocks
% PPD: 1x2 vector specifying  pixels per degree in x and y
%
% OUTPUT:
% events_eyelink
%
% Author: Nicolai Wolpert
% Email: Nicolai.Wolpert@ens.fr
% Version: 15.01.2018

fprintf('\n\n###############\nPreparing eyelink event structure...\n')


   
    events_eyelink = my_read_eyelink_asc(raw_data_path);

    % compute sample numbers for each event
    for i=1:length(events_eyelink.ssacc)

        % find sample number for current event
        [~, events_eyelink.ssacc(i).samplebegin] = min(abs((events_eyelink.dat(1,:)-events_eyelink.ssacc(i).timestampbegin)));

    end
    for i=1:length(events_eyelink.esacc)

        % find sample number for current event
        [~, events_eyelink.esacc(i).samplebegin] = min(abs((events_eyelink.dat(1,:)-events_eyelink.esacc(i).timestampbegin)));
        [~, events_eyelink.esacc(i).sampleend] = min(abs((events_eyelink.dat(1,:)-events_eyelink.esacc(i).timestampend)));

        % recalculate saccade amplitude in degrees based on the correct PPD
%         x_amplitude_pixels = events_eyelink.dat(2, events_eyelink.esacc(i).sampleend)-events_eyelink.dat(2, events_eyelink.esacc(i).samplebegin);
%         y_amplitude_pixels = events_eyelink.dat(3, events_eyelink.esacc(i).sampleend)-events_eyelink.dat(3, events_eyelink.esacc(i).samplebegin);
%         x_amplitude_degrees = x_amplitude_pixels/PPD(1);
%         y_amplitude_degrees = y_amplitude_pixels/PPD(2);
%         events_eyelink.esacc(i).saccamplitude = sqrt(x_amplitude_degrees.^2+y_amplitude_degrees.^2);
    end
    for i=1:length(events_eyelink.sblink)

        % find sample number for current event
        [~, events_eyelink.sblink(i).samplebegin] = min(abs((events_eyelink.dat(1,:)-events_eyelink.sblink(i).timestampbegin)));      

    end
    for i=1:length(events_eyelink.eblink)

        % find sample number for current event
        [~, events_eyelink.eblink(i).samplebegin] = min(abs((events_eyelink.dat(1,:)-events_eyelink.eblink(i).timestampbegin)));
        [~, events_eyelink.eblink(i).sampleend] = min(abs((events_eyelink.dat(1,:)-events_eyelink.eblink(i).timestampend)));

    end



end

