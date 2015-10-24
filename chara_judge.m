% This section of the code calculates which property from the previous
% calculation should be used in clustering process

% The assumption in this section is that the property number for one kind of
% wave should be the close. When we plot a histogram for one property, we
% should be able to get a rough idea of how many kinds of waves are there.
% To eliminate the problem of some unrelated properties influencing the
% performance of the clustering algorithm (for one channel), we tried to find which
% properties can help us to tell waveforms apart, and then use them to
% clustering. The properties are normalized to make sure that they are in
% the same scale when using k-means method.

% Manipulate data on its backup
totalcharadata = waveformproperty;
totalwavedata = reversedata;
originalwavedata = originaldata;
channelnumber = size(totalcharadata, 2);
for channel = 1:channelnumber
    % Find the max and min for each individual property inside each channel
    propertymaxvalue{channel} = max(totalcharadata{channel},[],1);
    propertyminvalue{channel} = min(totalcharadata{channel},[],1);
end
% Find the total number of properties (- 2 to be exact), as the first one and
% second one are not properties
propertynum = size(propertymaxvalue{1}, 2);
for channel = 1: channelnumber
    wavenum = size(totalcharadata{channel}, 1);
    for property = 3: propertynum
        % Normalize each properties to a score
        maxvalue = propertymaxvalue{channel}(1, property);
        minvalue = propertyminvalue{channel}(1, property);
        % The difference between max and min
        if abs(maxvalue) >= abs(minvalue)
            normalizeconst = abs(maxvalue);
            if (abs(maxvalue) == abs(minvalue)) && (abs(maxvalue) == 0)
                normalizeconst = 1;
            end
        else
            normalizeconst = abs(minvalue);
        end
        for wave = 1: wavenum
            % According to each wave's each property's distance from 0,
            % give them a score from -100 to 100
            distancematrix{channel}(wave,1) = totalcharadata{channel}(wave,1);
            distancematrix{channel}(wave,2) = totalcharadata{channel}(wave,2);
            distancematrix{channel}(wave,property) = (totalcharadata{channel}(wave,property) - 0) ./ normalizeconst .* 100;
            roundeddistancematrix{channel}(wave,property) = ceil(distancematrix{channel}(wave,property));
        end
        % Find the amount of waves that have a specific score (i.e. create
        % a histogram)
        countmatrix{channel}(property, 1) = 0;
        for i = -100: 100
            countmatrix{channel}(property, i + 102) = size(find(roundeddistancematrix{channel}(:, property) == i), 1);
        end
        countmatrix{channel}(property, 203) = 0;
        %
        % If there are too few waves with the specific score, just ignore
        % them
        %{
        newcountmatrix{channel}(property, :) = countmatrix{channel}(property, :)/wavenum;
        figure
        plot(1:103, newcountmatrix{channel}(property, :));
        %}
        %
        %{
        for i = 1:103
            if countmatrix{channel}(property,i) < 0.015 * wavenum
                countmatrix{channel}(property,i) = 0;
            end
        end
        %}
        
        % Plot out the histogram
        %{
        newcountmatrix{channel}(property, :) = countmatrix{channel}(property, :)/wavenum;
        figure
        plot(1:203, newcountmatrix{channel}(property, :));
        %}
        
        % This section is one previously tested idea.
        % The idea here is to match the distribution of properties into
        % several normal distributions. However, I didn't find a good
        % implemented algorithm in terms of matching into one distribution.
        % Here is a simple implementation of this idea, however it is naive
        % and didn't work so well.
        %{
        % Combine individual sections that are seperated from each
        % other by 0 waves with that specific score (The idea here is
        % similar to trying to match each score distribution peak into a normal 
        % distribution)
        % For example, if the histogram is 0 1 2 0 3, then 1 and 2 will be
        % grouped into one section
        peaknum = 0;
        location = 1;
        % Find the start and end position for each section
        while location <=202
            if countmatrix{channel}(property,location) == 0 && countmatrix{channel}(property,location + 1) ~= 0
                peaknum = peaknum + 1;
                start_location{peaknum} = location + 1;
            elseif countmatrix{channel}(property,location) ~= 0 && countmatrix{channel}(property,location + 1) == 0
                end_location{peaknum} = location;
            end
            location = location + 1;
        end
        % Sum each section's total number to a single point
        new_countmatrix(1, 1:203) = 0;
        for i = 1: peaknum
            for j = start_location{i} : end_location{i}
                new_countmatrix(1, ceil((start_location{i}+end_location{i})/2)) = countmatrix{channel}(property,j) + new_countmatrix(1, ceil((start_location{i}+end_location{i})/2));
            end
        end
        % countmatrix{channel}(property,:) = new_countmatrix(1, :);
        %}
    end
end

% With the distribution data, we can see which properties are significant
% and which can help us in clustering
for channel = 1: channelnumber
    wavenum = size(totalcharadata{channel}, 1);
    for property = 3:propertynum
        % Here is another idea. we look at how many peaks are there in the
        % histogram, so as to determine how many normal distributions are
        % there in the property data.
        
        % Find the local maximum for the distribution data for a single
        % property
        [maximavalue{channel}{property}, maxindices{channel}{property}] = findpeaks(countmatrix{channel}(property,:));
        if ~isempty(maximavalue{channel}{property})
            if size(maximavalue{channel}{property},2) >= 3
                % Calculate the second-order maximum to make sure that we
                % found the close-to-real maximum
                [twicemaximavalue{channel}{property}, twicemaxindices{channel}{property}] = findpeaks(maximavalue{channel}{property});
                largeenoughnumber{channel}{property} = size(twicemaximavalue{channel}{property},2);
            elseif size(maximavalue{channel}{property},2) == 2
                largeenoughnumber{channel}{property} = 2;
            else
                largeenoughnumber{channel}{property} = 1;
            end
        else
            largeenoughnumber{channel}{property} = 0;
        end
        %{
        % Find local maximums that are large enough for consisting a specific
        % pattern of spike
        largeenoughindices{channel}{property} = find(maximavalue{channel}{property} > (0.2 .* wavenum));
        % How many patterns can we find from this property
        largeenoughnumber{channel}{property} = size(largeenoughindices{channel}{property},2);
        %}
        %{
        % Find the distance between each significant "normal distribution's peak"
        differencenum = 1;
        if largeenoughnumber{channel}{property} > 1
            for i = 1: (largeenoughnumber{channel}{property} - 1)
                for j = 1: (largeenoughnumber{channel}{property} - i)
                    indiciesdifference{channel}{property}(differencenum,1) = abs(maxindices{channel}{property}(1,largeenoughindices{channel}{property}(1,i)) - maxindices{channel}{property}(1,largeenoughindices{channel}{property}(1,i + j)));
                    differencenum = differencenum + 1;
                end
            end
        else
            indiciesdifference{channel}{property}(1,1) = 0;
        end
        %}
        % There should be more than two peaks in the histogram so as to
        % consider the property being useful for clustering.
        if (largeenoughnumber{channel}{property} > 1)
            propertyuseornot(channel, property) = 1;
        else
            propertyuseornot(channel, property) = 0;
        end
    end
    % Property 1 and 2 will always be used.
    propertyuseornot(:, 1) = 1;
    propertyuseornot(:, 2) = 1;
end
% Grab the properties that are significant, and put them into a matrix that will
% be used for clustering
for channel = 1:channelnumber
    clusteringpropertymatrix{channel} = distancematrix{channel}(:, find(propertyuseornot(channel,:)));
end
% propertycomparematrix = totalcharadata;
propertycomparematrix = clusteringpropertymatrix;
%clusteringpropertymatrix = distancematrix;
%{
propertycomparematrix = distancematrix;
for channel = 1:channelnumber
    propertycomparematrix{channel}(:,1) = totalcharadata{channel}(:, 1);
    propertycomparematrix{channel}(:,2) = totalcharadata{channel}(:, 2);
end
%}