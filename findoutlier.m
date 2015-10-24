function [outlierdata, outlierproperty, outlierindex] = findoutlier(wavedata, propertydata)
% Find outliers inside each group by comparing properties and shape
pd = makedist('Normal');
multiplier = icdf(pd, 0.975); % Two sided 0.95 propability
numgroup = size(wavedata,2);
propertynum = size(propertydata{1},2) - 2;
propertynum;
outlierdata = [];
outlierproperty = [];
outlierindex = [];
% Build a 95% confidence interval
for i = 1: numgroup
    propertystd(i,:) = std(propertydata{i}(:, 3:(propertynum + 2)),1);
    propertymean(i,:) = mean(propertydata{i}(:, 3: (propertynum + 2)),1);
    for j = 1:propertynum
        limit{i}(1,j) = propertymean(i,j) + multiplier .* propertystd(i,j);
        limit{i}(2,j) = propertymean(i,j) - multiplier .* propertystd(i,j);
    end
    waveformmean{i} = calcmeanwaveform(wavedata{i}, propertydata{i});
end
% Calculate how many properties for one spike are not within the 95%
% confidence interval, as well as the waveform difference std
for i = 1: numgroup
    numingroup = size(propertydata{i},1);
    for j = 1: numingroup
        notequal{i}(j,1) = 0;
        for k = 1: propertynum
            if ((propertydata{i}(j, k + 2) >= limit{i}(1,k))||(propertydata{i}(j, k + 2) <= limit{i}(2,k)))
                notequal{i}(j,1) = notequal{i}(j,1) + 1;
            end
        end
        i;
        waveformdifference{i}(j,:) = wavedata{i}(j,:) - waveformmean{i}(propertydata{i}(j,1) + 1,:);
    end
    waveformdifferencestd{i} = std(waveformdifference{i},0,2);
end
clear relativeprop
for i = 1: numgroup
    notequalmean{i} = mean(notequal{i}, 1);
    notequalstd{i} = std(notequal{i}, 0, 1);
    % Make a Poisson distribution for the number of properties that differs
    poissonpd = makedist('Normal','mu',notequalmean{i},'sigma',notequalstd{i});
    poissonthreshold = icdf(poissonpd, 0.85);
    differencestdmean{i} = mean(waveformdifferencestd{i}, 1);
    differencestdstd{i} = std(waveformdifferencestd{i}, 0, 1);
    % Make a normal distribution for the waveshape differnece std
    normalpd = makedist('Normal','mu',differencestdmean{i},'sigma',differencestdstd{i});
    % Find those spikes whose significantly-different properties numbers
    % and waveshape difference std is large
    normalthreshold = icdf(normalpd, 0.85);
    propertyprobabilyoutlier{i} = find(notequal{i} >= poissonthreshold);
    waveformprobabilyoutlier{i} = find(waveformdifferencestd{i} >= normalthreshold);
    probabilyoutlier{i}=[];
    % Combine all the outliers into a single set
    while (size(propertyprobabilyoutlier{i},1) || size(waveformprobabilyoutlier{i},1))
        if ~size(propertyprobabilyoutlier{i},1)
            probabilyoutlier{i} = [probabilyoutlier{i};waveformprobabilyoutlier{i}(:,1)];
            waveformprobabilyoutlier{i}(:,1) = [];
            break
        elseif ~size(waveformprobabilyoutlier{i},1)
            probabilyoutlier{i} = [probabilyoutlier{i};propertyprobabilyoutlier{i}(:,1)];
            propertyprobabilyoutlier{i}(:,1) = [];
            break
        end
        if propertyprobabilyoutlier{i}(1,1) > waveformprobabilyoutlier{i}(1,1)
            probabilyoutlier{i} = [probabilyoutlier{i};waveformprobabilyoutlier{i}(1,1)];
            waveformprobabilyoutlier{i}(1,:) = [];
        elseif propertyprobabilyoutlier{i}(1,1) < waveformprobabilyoutlier{i}(1,1)
            probabilyoutlier{i} = [probabilyoutlier{i};propertyprobabilyoutlier{i}(1,1)];
            propertyprobabilyoutlier{i}(1,:) = [];
        else
            probabilyoutlier{i} = [probabilyoutlier{i};propertyprobabilyoutlier{i}(1,1)];
            propertyprobabilyoutlier{i}(1,:) = [];
            waveformprobabilyoutlier{i}(1,:) = [];
        end
    end
    clear relativeprop
    clear relativetotalprop
    outliernum{i} = size(probabilyoutlier{i},1);
    % Calculate the data that will be outputed for further use
    if outliernum{i}
        outlierdata = [outlierdata;wavedata{i}(probabilyoutlier{i}(:,1),:)];
        outlierproperty = [outlierproperty;propertydata{i}(probabilyoutlier{i}, :)];
        outlierindex{i} = probabilyoutlier{i};
    end
end