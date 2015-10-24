% With the property matrix that can be used in clustering, this part of the
% code is used to do the clustering
% However, for finding outliers, all the properties are used because right
% now those properties would have a smaller influence compared with
% the machine learning clustering process

clear X
clear idx
clear catgory
clearvars -except file noisermsmean originalaverage clusteringpropertymatrix propertycomparematrix totalwavedata originaldata reversedata waveformproperty reversepositivedata originalpositivedata positivewaveformproperty originalstd posoriginalstd sizewindow resultpropertymatrix
THRESHOLD = 1/((noisermsmean) - 0.025); % Trial and error result, 0.011
channels = size(originaldata,2);
for chn = 1:channels
    maxunits = 7;
    outlierdata{chn} = [];
    outlierproperty{chn} = [];
    outlieroriginaldata{chn} = [];
    round = 0;
    propertynum = size(clusteringpropertymatrix{chn},2);
    % If there are only two properties in our clustering matrix, it means
    % that it contains only the first two entry, so we will use the whole
    % property matrix to clustering
    if propertynum ==2
        clusteringpropertymatrix{chn} = propertycomparematrix{chn};
    end
    propertynum = size(clusteringpropertymatrix{chn},2);
    X = clusteringpropertymatrix{chn}(:,3:propertynum);
    numwaves = size(X,1);
    ro = 1:size(reversedata{chn},2);
    clear idx
    clear sumd
    % Clustering from two 2 clusters to 7 clusters 
    [idx{2},ctrs{chn}{2}, sumd] = kmeans(X,2,...
        'Distance','city',...
        'Replicates',10);
    sumsumd = sum(sumd);
    for numofcluster = 3:maxunits;
        [idx{numofcluster},ctrs{chn}{numofcluster}, sumd] = kmeans(X,numofcluster,...
        'Distance','city',...
        'Replicates',10);
        sizemax = maxunits;
    end
    putintocatgory = 1; % A boolean variable to tell if there still exists 
    % clusters that are too small
    while putintocatgory
        clear tempcatgory
        clear tempproperty
        clear tempclusteredproperty
        clear temporiginalwave
        clear temptimestampgroup
        clear indexes
        for j = 1: sizemax
            indexes(j)=1;
        end
        % Base on the clustering result, put each spikes into their
        % specific groups
        for i = 1: numwaves
            for j = 1:sizemax
                if (idx{sizemax}(i) == j)
                    tempcatgory{j}(indexes(j),:) = totalwavedata{chn}(i,:);
                    tempproperty{j}(indexes(j),:) = propertycomparematrix{chn}(i,:);
                    temporiginalwave{j}(indexes(j),:) = originaldata{chn}(i,:);
                    tempclusteredproperty{j}(indexes(j),:) = clusteringpropertymatrix{chn}(i,:);
                    indexes(j) = indexes(j) + 1;
                end
            end
        end
        % If there exists group with less than 30 spikes in it, use the
        % clustering result with fewer groups, and then repeat the while
        % loop
        minvalue = min(indexes);
        if minvalue < 30
            sizemax = sizemax - 1;
        elseif minvalue >= 30
            putintocatgory = 0;
        end
    end
    % Copy the clustering result into a variable that will be further
    % passed to other parts of the code
    originalwave{chn} = temporiginalwave;
    catgory{chn} = tempcatgory;
    property{chn} = tempproperty;
    clusteredproperty{chn} = tempclusteredproperty;
    inter_exist = 1; % Boolean variable, whether there exist groups that
    % are very similar and can be combines
    times = 0;
    while (inter_exist)
        numgroup = size(catgory{chn},2)
        propertynum = size(clusteredproperty{chn}{numgroup},2) - 2;
        % Call finding outlier function to find outliers in the group
        clear tempoutlieroriwave
        tempoutlieroriwave = [];
        if (round < 7)
            [tempoutlierdata, tempoutlierproperty, tempoutlierindex] = findoutlier(catgory{chn},property{chn});
            % Combine all the outliers from one channel into a set (group)
            for i = 1: numgroup
                tempoutlieroriwave = [tempoutlieroriwave;originalwave{chn}{i}(tempoutlierindex{i}(:,1), :)];
            end
            % Clear those outliers from the clustering result
            for i = 1:numgroup
                originalwave{chn}{i}(tempoutlierindex{i},:) = [];
                catgory{chn}{i}(tempoutlierindex{i},:) = [];
                property{chn}{i}(tempoutlierindex{i},:) = [];
                clusteredproperty{chn}{i}(tempoutlierindex{i},:) = [];
            end
            % Put outliers into a set (group) that contains previous steps'
            % outliers (all the outliers for this channel)
            outlierdata{chn} = [outlierdata{chn};tempoutlierdata];
            outlieroriginaldata{chn} = [outlieroriginaldata{chn};tempoutlieroriwave];
            outlierproperty{chn} = [outlierproperty{chn};tempoutlierproperty];
        end
        
        clear tempnewcatgory
        clear tempnewproperty
        clear tempnewclusteringproperty
        clear totalthrehsold
        clear tempneworiginalwave
        round = round + 1;
        ro = 1:size(reversedata{chn},2); % The x-axis's range that will be 
        % used for plotting graph
        
        % Compare the shape and properties of different groups using
        % statistics that are implemented in the calcstdratio function.
        % Offer values that will be used to determine if two groups are
        % very similar and can be combined
        [tempshapethreshold, tempwavediffgroupstdmean, tempwavediffgroupstdstd, tempwavedifferencestdmean, tempwavedifferencestdstd] = calcstdratio(catgory{chn}, property{chn});
        temppropertythreshold = calcpropertyratio(property{chn});
        % Scan through the whole value table to find the min value and its
        % position
        changed = 0; % Boolean variable
        mintotal = (tempshapethreshold(1,1) + temppropertythreshold(1,1)) ./ 2;
        mintotalrow = 1;
        mintotalcol = 1;
        for i = 1: (numgroup - 1)
            if(numgroup - i) >0
                for j = 1: (numgroup - i)
                    % Average the calculated value between shape and
                    % property, and use this value to determine similarity
                    totalthrehsold(i,j) = (tempshapethreshold(i,j) + temppropertythreshold(i,j)) ./ 2;
                  
                    if (totalthrehsold(i,j) < mintotal)
                        mintotal = totalthrehsold(i,j);
                        mintotalrow = i;
                        mintotalcol = j;
                    end
                end
            else
                break
            end
        end
        % Store the calculated value for debug purpose
        totalpermthrehsold{chn} = totalthrehsold;
        % If the min value from the table is smaller than the set
        % threshold, combine the two groups
        if (mintotal <= THRESHOLD)
            changed = 1;
        end
        if changed
            %If there are only two groups, do not combine them
            if numgroup == 2
                break
            end
        else
            % If not going to combine, do the same calculation for similarity, and
            % save the result for debug purpose
            [tempshapethreshold, tempwavediffgroupstdmean, tempwavediffgroupstdstd, tempwavedifferencestdmean, tempwavedifferencestdstd] = calcstdratio(catgory{chn}, property{chn});
            clear totalthrehsold
            temppropertythreshold = calcpropertyratio(property{chn});
            for i = 1: (numgroup - 1)
                if(numgroup - i) >0
                    for j = 1: (numgroup - i)
                        totalthrehsold(i,j) = (tempshapethreshold(i,j) + temppropertythreshold(i,j)) ./ 2;
                    end
                else
                    break
                end
            end
            totalpermthrehsold{chn} = totalthrehsold;
            shapethreshold{chn} = tempshapethreshold;
            propertythreshold{chn} = calcpropertyratio(property{chn});
            wavediffgroupstdmean{chn} = tempwavediffgroupstdmean;
            wavediffgroupstdstd{chn} = tempwavediffgroupstdstd;
            wavedifferencestdmean{chn} = tempwavedifferencestdmean;
            wavedifferencestdstd{chn} = tempwavedifferencestdstd;
            break
        end
        inter_exist = 0; % Boolean variable that will be used to 
        % determine if further combination is needed
        % Calculate the total amount of spikes in the data
        wavenum = 0;
        for i = 1:numgroup
            wavenum = size(catgory{chn}{i},1) + wavenum;
        end
        row = mintotalrow;
        col = mintotalcol;
        if (changed)
            % Based on the min value, combine groups
            for j = 1:numgroup
                if j < row + col
                    tempnewcatgory{j} = catgory{chn}{j};
                    tempnewproperty{j} = property{chn}{j};
                    tempneworiginalwave{j} = originalwave{chn}{j};
                    tempnewclusteringproperty{j} = clusteredproperty{chn}{j};
                elseif j == row + col
                    tempnewcatgory{row} = [tempnewcatgory{row};catgory{chn}{j}];
                    tempnewproperty{row} = [tempnewproperty{row};property{chn}{j}];
                    tempneworiginalwave{row} = [tempneworiginalwave{row};originalwave{chn}{j}];
                    tempnewclusteringproperty{row} = [tempnewclusteringproperty{row};clusteredproperty{chn}{j}];
                else
                    tempnewcatgory{j - 1} = catgory{chn}{j};
                    tempnewproperty{j - 1} = property{chn}{j};
                    tempneworiginalwave{j - 1} = originalwave{chn}{j};
                    tempnewclusteringproperty{j - 1} = clusteredproperty{chn}{j};
                end
            end
            
            %{
            % Calculate the std value before combination and after
            % combination. This step is used to see if we have dramatically
            % increased the std value (likely wrongly combined two groups
            % that are not similar at all
            stdstdold = calcstd(catgory{chn}, property{chn});
            stdstdnew = calcstd(tempnewcatgory, tempnewproperty);
            
            % If the ratio is larger than 4 (randomly determined)
            % Do not combine them, and terminate the combining program
            if (stdstdnew(row,1) ./ stdstdold(row,1) > 4)
                 tempnewcatgory = catgory{chn};
                 tempnewproperty = property{chn};
                 tempneworiginalwave = originalwave{chn};
                 tempnewclusteringproperty = clusteredproperty{chn};
                 break
            end
              %}
            
            % If not, decrease the number of group
            numgroup = numgroup - 1;
            opts = statset('Display','final');
            clear ctrs_new
            % Find the new point where the center of cluster is
            for i = 1: numgroup
                X = tempnewclusteringproperty{i}(:, 3:(propertynum + 2));
                [idx_new,ctrs_new(i,:), sumd_new] = kmeans(X,1,...
                    'Distance','city',...
                    'Replicates',1);
            end
            ctrs{chn}(:, numgroup + 1) = [];
            
        else
            % If there are no similar groups, just keep the clustering result
            % as they are right now
            tempnewcatgory = catgory{chn};
            tempnewproperty = property{chn};
            tempnewclusteringproperty = clusteredproperty{chn};
            tempneworiginalwave = originalwave{chn};
            ctrs_new = ctrs{chn}{numgroup};
        end
        % Put the combining (or did not combine) result into the matrix
        % that will be carried down further
        catgory{chn} = tempnewcatgory;
        property{chn} = tempnewproperty;
        originalwave{chn} = tempneworiginalwave;
        clusteredproperty{chn} = tempnewclusteringproperty;
        ctrs{chn}{numgroup} = ctrs_new;
        clear totalthrehsold
        % Do the same calculation again to see if there exist similar
        % group, if so, will execute the combination step again (this part
        % of code will run two times currrently, can be combined into one
        % if we need to improve efficiency
        [tempshapethreshold, tempwavediffgroupstdmean, tempwavediffgroupstdstd, tempwavedifferencestdmean, tempwavedifferencestdstd] = calcstdratio(catgory{chn}, property{chn});
        temppropertythreshold = calcpropertyratio(property{chn});
        mintotal = (tempshapethreshold(1,1) + temppropertythreshold(1,1)) ./ 2;
        mintotalrow = 1;
        mintotalcol = 1;
        for i = 1: (numgroup - 1)
            if(numgroup - i) >0
                for j = 1: (numgroup - i)
                    totalthrehsold(i,j) = (tempshapethreshold(i,j) + temppropertythreshold(i,j)) ./ 2;
                    if (totalthrehsold(i,j) < mintotal)
                        mintotal = totalthrehsold(i,j);
                        mintotalrow = i;
                        mintotalcol = j;
                    end
                end
            else
                break
            end
        end
        totalpermthrehsold{chn} = totalthrehsold;
        if (mintotal <= THRESHOLD)
            inter_exist = 1;
        end
        
        % Save for debug use
        shapethreshold{chn} = tempshapethreshold;
        propertythreshold{chn} = temppropertythreshold;
        wavediffgroupstdmean{chn} = tempwavediffgroupstdmean;
        wavediffgroupstdstd{chn} = tempwavediffgroupstdstd;
        wavedifferencestdmean{chn} = tempwavedifferencestdmean;
        wavedifferencestdstd{chn} = tempwavedifferencestdstd;
        %{
        % Plot the clustering result after each combination
        figure
        for j = 1: numgroup
            subplot(2,3,j)
            wavenum = size(catgory{chn}{j},1);
            for i = 1: wavenum
                plot(ro, catgory{chn}{j}(i,:));
                hold on;
            end
            %axis([0 100 -1.5 1.5])
        end
        a=1;
        %}
        a=1;
    end
    %{
    % Plot the final clustering result
    figure
    for j = 1: numgroup
        subplot(2,3,j)
        wavenum = size(catgory{chn}{j},1);
        for i = 1: wavenum
            plot(ro, catgory{chn}{j}(i,:));
            hold on;
        end
        a=1;
        axis([0 76 -1.5 1.5])
        %axis([0 50 -0.00005 0.00005])
    end
    %}
end