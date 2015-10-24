% In this step, outliers are compared to each group and put into a group
% that they may very likely belong to
% In this section, the assumption is that the amount of
% significantly-different properties should follow a normal distrbution
% compared to another group. Also the waveshape difference's standard
% deviation should also follow a normal distribution

% First, process the data to find statistical mean and std data for each
% properties inside each group
channels = size(catgory,2);
pd = makedist('Normal');
multiplier = icdf(pd, 0.975); % Two sided 0.95 propability
for chn =1:channels;
    chn;
    numgroup = size(catgory{chn},2);
    propertynum = size(property{chn}{1},2) - 2;
    for i = 1: numgroup
        propertystd(i,:) = std(property{chn}{i}(:, 3:(propertynum + 2)),1);
        propertymean(i,:) = mean(property{chn}{i}(:, 3: (propertynum + 2)),1);
        
        % Calculate the 95% confidence interval for each property
        for j = 1:propertynum
            limit{i}(1,j) = propertymean(i,j) + multiplier .* propertystd(i,j);
            limit{i}(2,j) = propertymean(i,j) - multiplier .* propertystd(i,j);
        end
        
        % Calculate the mean waveform for each cluster
        waveformmean{i} = calcmeanwaveform(catgory{chn}{i}, property{chn}{i});
    end
    for i = 1: numgroup
        numingroup = size(property{chn}{i},1);
        for j = 1: numingroup
            notequal{chn}{i}(j,1) = 0;
            for k = 1: propertynum
                % Calculate how many properties are not within the range of
                % the confidence interval
                if ((property{chn}{i}(j, k + 2) >= limit{i}(1,k))||(property{chn}{i}(j, k + 2) <= limit{i}(2,k)))
                    notequal{chn}{i}(j,1) = notequal{chn}{i}(j,1) + 1;
                end
            end
            % Calculate the difference between mean waveform and each
            % spikes
            waveformdifference{chn}{i}(j,:) = catgory{chn}{i}(j,:) - waveformmean{i}(property{chn}{i}(j,1) + 1,:);
        end
        % Calculate the standard deviation for the above difference
        waveformdifferencestd{chn}{i} = std(waveformdifference{chn}{i},0,2);
    end
    clear relativeprop
    %figure
    for i = 1: numgroup
        % Arrange existing data into matrixes that will be used in the
        % rearranging outlier part
        notequalmean{chn}{i} = mean(notequal{chn}{i}, 1);
        notequalstd{chn}{i} = std(notequal{chn}{i}, 0, 1);
        originalcatgory{chn}{i} = catgory{chn}{i};
        originalproperty{chn}{i} = property{chn}{i};
        originalnotequal{chn}{i} = notequal{chn}{i};
        originalwaveshape{chn}{i} = originalwave{chn}{i};
        originalwaveformdifferencestd{chn}{i} = waveformdifferencestd{chn}{i};
        a = size(originalcatgory{chn}{i},1);
        %{
        for k = 1: a
            subplot(2,2,i)
            plot(ro, originalcatgory{chn}{i}(k,:))
            hold on
            axis([0 100 -1.5 1.5])  
        end
        %}
    end
    clear propertystd propertymean
end

% Then we will start to rearrange outliers
channels = size(outlierproperty,2);
pd = makedist('Normal');
multiplier = icdf(pd, 0.975); % Two sided 0.95 propability
relativeprop = [];
propertynum = size(originalproperty{1}{1},2);
clear poissonpd
clear normalpd
for chn = 1: channels
    groupnum = size(originalproperty{chn},2);
    propertynum = size(originalproperty{chn}{1},2);
    % For each cluster, calculate a 95% confidence interval for each of its
    % properties, as well as the mean waveform difference standard
    % deviation
    for group = 1: groupnum
        origipropertymean{chn}{group} = mean(originalproperty{chn}{group},1);
        originalwaveformmean{chn}{group} = calcmeanwaveform(originalcatgory{chn}{group},originalproperty{chn}{group});
        origipropertystd{chn}{group} = std(originalproperty{chn}{group},1);
        for property = 1:propertynum
            boundary{chn}{group}(1,property) = origipropertymean{chn}{group}(1,property) + multiplier .* origipropertystd{chn}{group}(1,property);
            boundary{chn}{group}(2,property) = origipropertymean{chn}{group}(1,property) - multiplier .* origipropertystd{chn}{group}(1,property);
        end
        % What is the mean number and its std of properties that are different
        % from the cluster mean
        originalnotequalmean{chn}{group} = mean(originalnotequal{chn}{group},1);
        originalnotequalstd{chn}{group} = std(originalnotequal{chn}{group},1);
        % What is the mean value and std for the waveform shape difference
        % std for each individual spikes shape compared with the mean shape
        originalwaveformdifferencestdmean{chn}{group} = mean(originalwaveformdifferencestd{chn}{group},1);
        originalwaveformdifferencestdstd{chn}{group} = std(originalwaveformdifferencestd{chn}{group},1);
        %poissonpd{chn}{group} = makedist('Poisson','lambda',originalnotequalmean{chn}{group});
        % Make the number of not-equal properties a normal distribution
        poissonpd{chn}{group} = makedist('Normal','mu',originalnotequalmean{chn}{group},'sigma',originalnotequalstd{chn}{group});
        % Make the shape difference std a normal distribution
        normalpd{chn}{group} = makedist('Normal','mu',originalwaveformdifferencestdmean{chn}{group},'sigma',originalwaveformdifferencestdstd{chn}{group});
    end
end
newcatgory = originalcatgory;
newproperty = originalproperty;
neworiginalshape = originalwaveshape;
for chn = 1: channels
    chn;
    outlierdelete{chn} = [];
    outliernum{chn} = size(outlierproperty{chn},1);
    groupnum = size(originalproperty{chn},2);
    propertynum = size(originalproperty{chn}{1},2);
    % For each outlier, compare it to existing groups, calculate the
    % number of properties that do not belong to the group, as well as the
    % shape difference std
    for i = 1: outliernum{chn}
        for j = 1: groupnum
            % Calculate the probability
            relativelesslikelynum{chn}{i}(j,1) = 0;
            relativelesslikelynum{chn}{i}(j,2) = j;
            differencestd{chn}{i}(j,2) = j;
            for m = 3: propertynum
                if (outlierproperty{chn}(i,m) >= boundary{chn}{j}(1,m)) || (outlierproperty{chn}(i,m) <= boundary{chn}{j}(2,m))
                    relativelesslikelynum{chn}{i}(j,1) = relativelesslikelynum{chn}{i}(j,1) + 1;
                end
            end
            waveformrelativedifferece{chn}{i}(j,:) = outlierdata{chn}(i,:) - originalwaveformmean{chn}{j}(outlierproperty{chn}(i,1) + 1,:);
            waveformrelativedifferecestd{chn}{i} = std(waveformrelativedifferece{chn}{i},0,2);
            % Calculate the normalized value for the two calculated number being
            % put into the cluster distribution
            numberprop{chn}{i}(j,1) = (relativelesslikelynum{chn}{i}(j,1) - originalnotequalmean{chn}{j}) / originalnotequalstd{chn}{j};
            diffprop{chn}{i}(j,1) = (waveformrelativedifferecestd{chn}{i}(j,1) - originalwaveformdifferencestdmean{chn}{j}) / originalwaveformdifferencestdstd{chn}{j};
            totalprop{chn}{i}(j,1) = j;
            % Average the two values
            totalprop{chn}{i}(j,2) = (numberprop{chn}{i}(j,1) + diffprop{chn}{i}(j,1))/2;
        end
        % Find the largest probability
        sortednum{chn}{i} = sortrows(totalprop{chn}{i},2);
        % If the smallest normalized value is smaller than 10, then we will
        % regard the outlier as belong to this group
        if sortednum{chn}{i}(1, 2) < 10
            newcatgory{chn}{sortednum{chn}{i}(1,1)} = [newcatgory{chn}{sortednum{chn}{i}(1,1)};outlierdata{chn}(i,:)];
            newproperty{chn}{sortednum{chn}{i}(1,1)} = [newproperty{chn}{sortednum{chn}{i}(1,1)};outlierproperty{chn}(i,:)];
            neworiginalshape{chn}{sortednum{chn}{i}(1,1)} = [neworiginalshape{chn}{sortednum{chn}{i}(1,1)};outlieroriginaldata{chn}(i,:)];;
            outlierdelete{chn} = [outlierdelete{chn};i,sortednum{chn}{i}(1,:)];
        end
    end
    % Delete combined outlier from the outlier group
    if ~isempty(outlierdelete{chn})
        outlierproperty{chn}(outlierdelete{chn}(:,1),:) = [];
        outlierdata{chn}(outlierdelete{chn}(:,1),:) = [];
    end
end