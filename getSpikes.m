function [ pks ] = getSpikes(signal)
%%% Output the indices of a signal where there are peaks %%%
[b a] = butter(4,.05,'low'); %setting the lowpass at .1*Fs=33Hz
signal_trend= filtfilt(b,a,signal);

%Take the difference of the raw data from the low-passed data
%should select for high frequency activity
differential = abs(signal- signal_trend);

meanline = zeros(1,length(differential)) + mean(differential);
sigmas = meanline + 4*std(differential);

differential_suppressed = differential;
differential_suppressed(differential<meanline + 4*std(differential)) = 0;
[~,pks] = findpeaks(differential_suppressed);

%Remove any peaks that are less than the low-passed trend
pks(signal(pks)<signal_trend(pks))=[];


end

