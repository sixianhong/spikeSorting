function stdoutput = calcstd (catgory,property)
% Calculate the waveshape difference standard deviation's standard deviation bewteen one group 
% and all the other groups. This one is used to determine if we wrongly
% combined two dissimilar group
numgroup = size(catgory,2);
for i = 1: numgroup
	tempelementnum{i} = size(property{i},1);
	tempmeanelement(i, :) = mean(property{i},1);
	tempstdelement(i, :) = std(property{i},1);
	tempmeanwaveform{i} = calcmeanwaveform(catgory{i},property{i});
end
for group = 1: numgroup
	wavenumber = size(catgory{group},1);
	for wave = 1: wavenumber
        tempwavedifference{group}(wave,:) = catgory{group}(wave,:) - tempmeanwaveform{group}(property{group}(wave,1) + 1,:);
    end
	tempwavedifferencestd{group} = std(tempwavedifference{group}, 0 , 2);
	tempwavedifferencestdmean(group,1) = mean(tempwavedifferencestd{group},1);
	tempwavedifferencestdstd(group,1) = std(tempwavedifferencestd{group}, 0 ,1);
end
stdoutput = tempwavedifferencestdstd;