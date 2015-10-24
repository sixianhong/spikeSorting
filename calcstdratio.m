function [threshold, tempwavediffgroupstdmean, tempwavediffgroupstdstd, tempwavedifferencestdmean, tempwavedifferencestdstd] = calcstdratio (catgory,property)
% Calculate the normalized score of waveshape difference standard deviation result bewteen one group 
% and all the other groups
% Input: waveform and property matrix
% Output: Normalized score for comparing two different groups' wave shape
numgroup = size(catgory,2);
% Calculate the mean waveform for each group
for i = 1: numgroup
	tempelementnum{i} = size(property{i},1);
	tempmeanelement(i, :) = mean(property{i},1);
	tempstdelement(i, :) = std(property{i},1);
	tempmeanwaveform{i} = calcmeanwaveform(catgory{i},property{i});
end
% Compare the group's own spikes's waveform to group mean shape, calculate 
% the standard deviation for the difference between spike data and mean waveshape 
for group = 1: numgroup
	wavenumber = size(catgory{group},1);
	for wave = 1: wavenumber
        tempwavedifference{group}(wave,:) = catgory{group}(wave,:) - tempmeanwaveform{group}(property{group}(wave,1) + 1,:);
    end
	tempwavedifferencestd{group} = std(tempwavedifference{group}, 0 , 2);
	tempwavedifferencestdmean(group,1) = mean(tempwavedifferencestd{group},1);
	tempwavedifferencestdstd(group,1) = std(tempwavedifferencestd{group}, 0 ,1);
end
% Use other groups' spikes' data, compare to the mean of one group, calculate the waveshape difference std
for i = 1: (numgroup - 1)
	wavenumber = size(catgory{i},1);
	for j = 1: (numgroup - i)
        for wave = 1:wavenumber
            tempwavediffgroup{i}{j}(wave,:) = catgory{i}(wave,:) - tempmeanwaveform{i + j}(property{i}(wave,1) + 1,:);
        end
        tempwavediffgroupstd{i}{j} = std(tempwavediffgroup{i}{j},0,2);
        tempwavediffgroupstdmean{i}{j} = mean(tempwavediffgroupstd{i}{j},1);
        tempwavediffgroupstdstd{i}{j} = std(tempwavediffgroupstd{i}{j}, 0, 1);
    end
end
% Normalize other groups' difference std using the group's own spikes' statistic results
for i = 1: (numgroup - 1)
	for j = 1: (numgroup - i)
    	tempthreshold(i,j) = abs(tempwavedifferencestdmean(i,1) - tempwavediffgroupstdmean{i}{j}) / tempwavedifferencestdstd(i,1);
    end
end
threshold = tempthreshold;