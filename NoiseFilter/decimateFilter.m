% DECIMATE NOISE FILTER

% COPYRIGHT: 
%
%   MIT License
%   © 2023 Dario Chiantello <pipesoundlibrary@gmail.com>
%

% INPUTS:
% 
%   firFilterOrder: desired filter order
%   windowFilterMask: filter mask to decimate
%

% OUTPUTS:
%
%   firFilterOrder: decimated filter order
%   decimatedFilterPulseResponse: impulse response of the decimated filter
%

function [firFilterOrder, decimatedFilterPulseResponse] = decimateFilter(firFilterOrder, windowFilterMask)

    numberOfSamplesPerWindow = 2*length(windowFilterMask)-2;

    %check filter order
    %the number of decimating points must divide the window in an integer number of intervals
    %the number of intervals  must be even to include the nyquist frequency
    while mod(numberOfSamplesPerWindow, firFilterOrder + 1) ~=0 || mod(firFilterOrder + 1,2)~=0

        firFilterOrder = firFilterOrder+1;

    end
    
    %create FIR filter with fir2
    %filterPulseResponse = fir2(firFilterOrder+1,maskNormalisedFrequencyBase,windowFilterMask);

    %create FIR filter by mask decimation and antitransform
    numberOfDecimatingPoints        = floor((firFilterOrder+1)/2) + 1; %firFilterOrder+1 always even  
    decimatedMask                   = windowFilterMask(1:floor(length(windowFilterMask)/(numberOfDecimatingPoints-1)):end);
    decimatedFiterSpectrum          = [decimatedMask' decimatedMask(end-1:-1:2)']';
    decimatedFilterPulseResponse    = real(ifft(decimatedFiterSpectrum));  %it's already real!
    decimatedFilterPulseResponse    = circshift(decimatedFilterPulseResponse, numberOfDecimatingPoints-1);
    decimatedFilterPulseResponse    = decimatedFilterPulseResponse .* hamming(firFilterOrder+1);
    
    