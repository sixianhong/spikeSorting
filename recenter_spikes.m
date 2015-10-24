% Center the spikes to their max position or min position if there is a greater
% than 1 standard deviation defference between the peak amplitude and the
% trough amplitude

% This code is not used in the final version of the algorithm because the
% performance is not quite as expected.


channels = size(reversedata,2);

newwindow = 0.002;
halfnumofpoints = newwindow * samprate(1) / 2;
posrecenter = 0;
negrecenter = 0;
% find the difference between max and mean, and min and mean
for channel = 1: channels
    meandata{channel} = mean(reversedata{channel},2);
    [maxdata{channel}, maxindex{channel}] = max(reversedata{channel},[] ,2);
    [mindata{channel}, minindex{channel}] = min(reversedata{channel},[] ,2);
    maxdiff{channel} = maxdata{channel} - meandata{channel};
    mindiff{channel} = meandata{channel} - mindata{channel};
end

% recenter spikes to their max or min position
for channel = 1: channels
    waves = size(reversedata{channel}, 1);
    halfpoints = size(reversedata{channel}, 2) / 2;
    for wave = 1: waves
        if maxdiff{channel}(wave,1) > 0.6 * (maxdiff{channel}(wave,1) + mindiff{channel}(wave,1))
            newreversedata{channel}(wave,:) = reversedata{channel}(wave, (maxindex{channel}(wave,1) - halfnumofpoints):(maxindex{channel}(wave,1) + halfnumofpoints));
            posrecenter = posrecenter + 1;
        else
            newreversedata{channel}(wave,:) = reversedata{channel}(wave, (minindex{channel}(wave,1) - halfnumofpoints): (minindex{channel}(wave,1) + halfnumofpoints));
            negrecenter = negrecenter + 1;
        end
    end
end

reversedata = newreversedata;