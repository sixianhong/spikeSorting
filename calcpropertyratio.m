function threshold = calcpropertyratio(propertydata)
% Calculate the normalized score of the number of properties that differ significantly between
% one group with all the other groups
% Input: property data matrix
% Output: Based on the calculation of the amount of properties that differ
% significantly from other group, calculate a normalized score that shows
% how different is the group from another group
groupnum = size(propertydata,2);
propertynum = size(propertydata{1},2);
pd = makedist('Normal');
multiplier = icdf(pd, 0.975);
for group = 1: groupnum
    propertymean{group} = mean(propertydata{group},1);
    propertystd{group} = std(propertydata{group},0,1);
    % Calculate a 95% confidence interval for each properties
    for property = 1: propertynum
        boundary{group}(1,property) = propertymean{group}(1,property) + multiplier .* propertystd{group}(1,property);
        boundary{group}(2,property) = propertymean{group}(1,property) - multiplier .* propertystd{group}(1,property);
    end
end
% For the target group, use its own spikes' properties
% Find the number of properties that is far away from the group mean
% for each spikes, then calculate a mean and std for the group
for group = 1:groupnum
    wavenum = size(propertydata{group},1);
    for wave = 1:wavenum
        numofproperty{group}(wave,1) = 0;
        for property = 3:propertynum
            % If the property is not within the confidence interval, treat
            % the property as an outlier, increase the number of properties
            % that differ significantly from the group mean
            if (propertydata{group}(wave,property) >= boundary{group}(1,property)) || (propertydata{group}(wave,property) <= boundary{group}(2,property))
                numofproperty{group}(wave,1) = numofproperty{group}(wave,1) + 1;
            end
        end
    end
    % Calculate the mean and std of the number of properties that differ
    % significantly from the group mean
    numofpropertymean{group} = mean(numofproperty{group},1);
    numofpropertystd{group} = std(numofproperty{group},0,1);
end
% For other groups, use their spikes' properties, compare to the target
% group's property confidence interval, and find property "outlier" number
for group = 1:(groupnum - 1)
    wavenum = size(propertydata{group},1);
    for tobecomparedgroup = 1: (groupnum - group) 
        for wave = 1:wavenum
            numofpropertycomparetoothergroup{group}{tobecomparedgroup}(wave,1) = 0;
            for property = 3:propertynum
                if (propertydata{group}(wave,property) >= boundary{group + tobecomparedgroup}(1,property)) || (propertydata{group}(wave,property) <= boundary{group + tobecomparedgroup}(2,property))
                    % Use the properties from other group's spikes to compare to the
                    % property confidence interval for each group, and see
                    % how many properties differ significantly (Similar to
                    % comparing the mean of two different groups but what
                    % we would get from here is the number of properties
                    % that differ significantly, rather than probability)
                    numofpropertycomparetoothergroup{group}{tobecomparedgroup}(wave,1) = numofpropertycomparetoothergroup{group}{tobecomparedgroup}(wave,1) + 1;
                end
            end
        end
        % Calculate a mean value for a group's all the spikes' number of properties that differ
        % significantly from the other group mean
        numofpropertycomparemean{group}{tobecomparedgroup} = mean(numofpropertycomparetoothergroup{group}{tobecomparedgroup},1);
    end
end
for group = 1:(groupnum - 1)
    for tobecomparedgroup = 1:(groupnum - group)
        % Normalize the number of different properties by the target group's mean
        % and std, and then we would get a normalized score
        threshold(group,tobecomparedgroup) = abs(numofpropertycomparemean{group}{tobecomparedgroup} - numofpropertymean{group}) ./ numofpropertystd{group};
    end
end