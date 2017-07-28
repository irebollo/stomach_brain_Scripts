function [filter_result]=tools_bpFilter(data,srate,filterOrder,center_frequency,filter_frequency_spread,transition_width,filterType)

%{

% Forward and reverse filter the signal using filtfilt. This corrects for phase distortion introduced by a one-pass filter.
% filterType has to be specified as string input e.g. ’fir2’
% filterOrder: length of filter kernel in samples (total length =forder +1)
    % Lenght in samples must be smaller than 1/3  of the data to be filtered 
% center_frequency: in hz
% filter_frequency_spread: in hz
% transition_width: in normalized units 0 to 1
% Data should be timeXvoxels 


IR commented 27/06/2017

%}
nyquist= srate*0.5;  % Nyquist frequency

% ffrequencies:  six numbers indicating the gain of the filter at different
% parts of the filter kernel 
% 0, start of lower transition zone, 
% lower bound of bandpass, upperband of bandpassm nyquist frequency and
% scale them to the nyquist frequency so nyq = one, %Frequencies are expressed in % of nyquist

ffrequencies   = [ 0 (1-transition_width)*(center_frequency-filter_frequency_spread)...
    (center_frequency-filter_frequency_spread) (center_frequency+filter_frequency_spread)...
    (1+transition_width)*(center_frequency+filter_frequency_spread) nyquist ]/nyquist;
idealresponse  = [ 0 0 1 1 0 0 ];
filterweights  = eval([filterType '(filterOrder,ffrequencies,idealresponse)']); % obtain filter kernel
filter_result = filtfilt(filterweights,1,data); % apply filter

end


