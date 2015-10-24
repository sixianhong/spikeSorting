function meanwaveform = calcmeanwaveform(wavedata, propertydata)
% Calculate the mean waveform using the wave data
% The output would contain two waveform, the first one is 
% the first one with the max at the midpoint
% the second one with the min at the midpoint
%
% Previously because centering seems to be a problem, I implemented the
% algorithm of calculating mean waveform (two of them) one centers at the 
% min and the other one centers at the max
%
% Input: shape and property matrix
% Output: two mean waveforms, one centered at the max, and the other one
% centered at min

mid_at_min = find(propertydata(:,1) == 0);
% Calculate a mean of the whole signal, in case the signal needs to be
% extended when recentering the mean waveshape
meanvalue = mean(mean(wavedata,2),1);
mid_at_min_num = size(mid_at_min,1);
% Calculate the mean shape for centering at min
midpoint_min_mean_waveform = mean(wavedata(mid_at_min,:));
% Calculate the mean shape for centering at max
mid_at_max = find(propertydata(:,1) == 1);
mid_at_max_num = size(mid_at_max,1);
midpoint_max_mean_waveform = mean(wavedata(mid_at_max,:));
% Find the min position for the centering-at-max's mean waveform
[value,min_location_for_midpoint_max] = min(midpoint_max_mean_waveform);
% Solve the case when no spike centers at the max
if mid_at_max_num
    [value,midlocation] = max(midpoint_max_mean_waveform);
else
    [value,midlocation] = min(midpoint_min_mean_waveform);
end
% 
total_datapoint_num = size(midpoint_max_mean_waveform,2);
% Elongate the data to fit the expected length
if min_location_for_midpoint_max <= midlocation
    shift_num = midlocation - min_location_for_midpoint_max;
    shifted_wave(1, 1:total_datapoint_num) = meanvalue;
    shifted_wave(1, shift_num + 1: total_datapoint_num) = midpoint_max_mean_waveform(1: total_datapoint_num - shift_num);
else
    shift_num = min_location_for_midpoint_max - midlocation;
    shifted_wave(1, 1:total_datapoint_num) = meanvalue;
    shifted_wave(1, 1: (total_datapoint_num - shift_num)) = midpoint_max_mean_waveform(shift_num + 1: total_datapoint_num);
end
% Calculate the new mean waveform that centers at min
if mid_at_min_num && mid_at_max_num
    new_meanwaveform = (midpoint_min_mean_waveform .* mid_at_min_num + shifted_wave .* mid_at_max_num) ./(mid_at_max_num + mid_at_min_num);
elseif mid_at_min_num && (~mid_at_max_num)
    new_meanwaveform = midpoint_min_mean_waveform;
elseif (~mid_at_min_num) && (mid_at_max_num)
    new_meanwaveform = shifted_wave;
else
    new_meanwaveform(1,1:total_datapoint_num) = meanvalue;
end
meanwaveform(1,:) = new_meanwaveform;
% Shift the calculated mean shape to centered-at-max case
[value,max_location] = max(new_meanwaveform);
if max_location <= midlocation
    shift_num = midlocation - max_location;
    output_shifted_wave(1, 1:total_datapoint_num) = meanvalue;
    output_shifted_wave(1, shift_num + 1: total_datapoint_num) = new_meanwaveform(1: total_datapoint_num - shift_num);
else
    shift_num = max_location - midlocation;
    output_shifted_wave(1, 1:total_datapoint_num) = meanvalue;
    output_shifted_wave(1, 1: (total_datapoint_num - shift_num)) = new_meanwaveform(shift_num + 1: total_datapoint_num);
end
meanwaveform(2,:) = output_shifted_wave;