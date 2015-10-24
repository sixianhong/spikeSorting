% Calculate properties from the de-noised signal

Fs = 50000;                    % Sampling frequency, used to calculate fast
% Fourier transform waveform
T = 1/Fs;                     % Sample time
L = size(reversedata{1}, 2);                     % Length of signal
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
f = Fs/2*linspace(0,1,NFFT/2+1);
multiplier = 2;
%tic

clear largest
clear secondlargest
channelnumber = size(reversedata, 2);
falltimeindex = [];
peakthentrough = [];
for channel = 1:channelnumber
    % Calculate general properties like max/min/threshold for each of the 
    % channels to improve performance
    wavenumber = size(reversedata{(channel)},1);
    datapointnumber = size(reversedata{(channel)},2);
    [minvalue, minindex] = min(reversedata{(channel)},[],2);
    [maxvalue, maxindex] = max(reversedata{(channel)},[],2);
    meanvalue = mean(reversedata{(channel)},2);
    maxmagnitude = maxvalue-minvalue;
    onetenthmagnitude = 0.1 .* maxmagnitude;
    upperthreshold = maxvalue - onetenthmagnitude;
    lowerthreshold = minvalue + onetenthmagnitude;
    standarddeviation = std(reversedata{channel},0,2);
    troughthreshold = meanvalue - multiplier .* standarddeviation;
    peakthreshold = meanvalue + multiplier .* standarddeviation;
    start_point = 0;
    end_point = 0;
    for wave = 1: wavenumber
        % Find the previous peak and following peak of the min trough
        % through walking around the trough
        prevfind = 0;
        prevlocation = minindex(wave,1);
        endfind = 0;
        endlocation = minindex(wave,1);
        while ~prevfind
            if prevlocation == 1
                break;
            end
            if (reversedata{channel}(wave,prevlocation - 1) > reversedata{channel}(wave,prevlocation))
                prevlocation = prevlocation - 1;
            else
                prevfind = 1;
            end
        end
        minaround(1,1) = prevlocation;
        while ~endfind
            if endlocation == datapointnumber
                break
            end
            if (reversedata{channel}(wave,endlocation) < reversedata{channel}(wave,endlocation + 1))
                endlocation = endlocation + 1;
            else
                endfind = 1;
            end
        end
        minaround(2,1) = endlocation;
        % Find the leading and following trough for the max peak
        prevfind = 0;
        prevlocation = maxindex(wave,1);
        endfind = 0;
        endlocation = maxindex(wave,1);
        while ~prevfind
            if prevlocation == 1
                break;
            end
            if (reversedata{channel}(wave,prevlocation - 1) < reversedata{channel}(wave,prevlocation))
                prevlocation = prevlocation - 1;
            else
                prevfind = 1;
            end
        end
        maxaround(1,1) = prevlocation;
        while ~endfind
            if endlocation == datapointnumber
                break
            end
            if (reversedata{channel}(wave,endlocation) > reversedata{channel}(wave,endlocation + 1))
                endlocation = endlocation + 1;
            else
                endfind = 1;
            end
        end
        maxaround(2,1) = endlocation;
        
        % Determine whether to use the peak or trough to calculate the rise 
        % time and fall time by comparing their height and chosse the one
        % that is higher or deeper (choose the dominate local max/min)
        % In other words, the rise/fall time is calculated from one
        % rise/fall, rather than using the max/min of the whole wave
        
        if (maxvalue(wave,1) - originalaverage{channel} >= abs(originalaverage{channel} - minvalue(wave,1)))
            if (reversedata{channel}(wave,maxaround(1,1)) <= reversedata{channel}(wave,maxaround(2,1)))
                minindex(wave,1) = maxaround(1,1);
                minvalue(wave,1) = reversedata{channel}(wave,maxaround(1,1));
            else
                minindex(wave,1) = maxaround(2,1);
                minvalue(wave,1) = reversedata{channel}(wave,maxaround(2,1));
            end
        else
            if (reversedata{channel}(wave,minaround(1,1)) >= reversedata{channel}(wave,minaround(2,1)))
                maxindex(wave,1) = minaround(1,1);
                maxvalue(wave,1) = reversedata{channel}(wave,minaround(1,1));
            else
                maxindex(wave,1) = minaround(2,1);
                maxvalue(wave,1) = reversedata{channel}(wave,minaround(2,1));
            end
        end
        %% calculate rise/ fall time by walking through the region
        
        % Case when first min point, and then max point
        if (maxindex(wave,1) > minindex(wave, 1))
            peakthentrough{(channel)}(wave,1) = 1;
            height = maxvalue(wave,1) - minvalue(wave,1);
            % Find the rise time using 60% level of change
            lowerthreshold = minvalue(wave, 1) + 0.2 .* height;
            higherthreshold = maxvalue(wave, 1) - 0.2 .* height;
            for i = minindex(wave, 1) : maxindex(wave, 1)
                if (reversedata{channel}(wave, i) <= lowerthreshold) && (reversedata{channel}(wave, i + 1) > lowerthreshold)
                    startpoint = i + (lowerthreshold - reversedata{channel}(wave, i)) ./ (reversedata{channel}(wave, i + 1) - reversedata{channel}(wave, i));
                    % Linear approximation
                elseif (reversedata{channel}(wave, i) <= higherthreshold) && (reversedata{channel}(wave, i + 1) > higherthreshold)
                    endpoint = i + (higherthreshold - reversedata{channel}(wave, i)) ./ (reversedata{channel}(wave, i + 1) - reversedata{channel}(wave, i));
                end
            end
            risetime = endpoint - startpoint;
            
            prevmaxnotfound = 1;
            prevmaxloc = minindex(wave,1);
            % Find the previous local maximum by walking through
            while prevmaxnotfound
                if prevmaxloc > 1
                    if reversedata{channel}(wave, prevmaxloc) < reversedata{channel}(wave, prevmaxloc - 1)
                        prevmaxloc = prevmaxloc - 1;
                    else
                        prevmaxnotfound = 0;
                    end
                else
                    prevmaxloc = 1;
                    prevmaxnotfound = 0;
                end
            end
            % Calculate the fall time before the rising edge
            % The method here is the same as the one used to calculate rise
            % time
            lowerthreshold = minvalue(wave, 1) + 0.2 * (reversedata{channel}(wave, prevmaxloc) - minvalue(wave, 1));
            higherthreshold = reversedata{channel}(wave, prevmaxloc) - 0.2 * (reversedata{channel}(wave, prevmaxloc) - minvalue(wave, 1));
            for i = prevmaxloc : minindex(wave, 1)
                if (reversedata{channel}(wave, i) >= lowerthreshold) && (reversedata{channel}(wave, i + 1) < lowerthreshold)
                    endpoint = i + (lowerthreshold - reversedata{channel}(wave, i)) ./ (reversedata{channel}(wave, i + 1) - reversedata{channel}(wave, i));
                elseif (reversedata{channel}(wave, i) >= higherthreshold) && (reversedata{channel}(wave, i + 1) < higherthreshold)
                    startpoint = i + (higherthreshold - reversedata{channel}(wave, i)) ./ (reversedata{channel}(wave, i + 1) - reversedata{channel}(wave, i));
                end
            end
            firstfalltime = endpoint - startpoint;
            % Find the following local min point
            follminnotfound = 1;
            follminloc = maxindex(wave,1);
            while follminnotfound
                if follminloc < datapointnumber
                    if reversedata{channel}(wave, follminloc) > reversedata{channel}(wave, follminloc + 1)
                        follminloc = follminloc + 1;
                    else
                        follminnotfound = 0;
                    end
                else
                    follminloc = datapointnumber;
                    follminnotfound = 0;
                end
            end
            higherthreshold = maxvalue(wave, 1) - 0.2 * (maxvalue(wave, 1) - reversedata{channel}(wave, follminloc));
            lowerthreshold = maxvalue(wave, 1) - 0.8 * (maxvalue(wave, 1) - reversedata{channel}(wave, follminloc));
            % Calculate the fall time following the rising edge
            for i = maxindex(wave, 1) : follminloc
                if i == follminloc
                    endpoint = follminloc;
                    break
                end
                if (reversedata{channel}(wave, i) >= lowerthreshold) && (reversedata{channel}(wave, i + 1) < lowerthreshold)
                    endpoint = i + (lowerthreshold - reversedata{channel}(wave, i)) ./ (reversedata{channel}(wave, i + 1) - reversedata{channel}(wave, i));
                elseif (reversedata{channel}(wave, i) >= higherthreshold) && (reversedata{channel}(wave, i + 1) < higherthreshold)
                    startpoint = i + (higherthreshold - reversedata{channel}(wave, i)) ./ (reversedata{channel}(wave, i + 1) - reversedata{channel}(wave, i));
                end
            end
            secondfalltime = endpoint - startpoint;
            % Fall time is negative while rise time is positive
            risefalltime{channel}(wave,1) = 0 - firstfalltime;
            risefalltime{channel}(wave,2) = risetime;
            risefalltime{channel}(wave,3) = 0 - secondfalltime;
        else
            % Do the same process for a first max, followed by a min case
            peakthentrough{(channel)}(wave,1) = 0;
            height = maxvalue(wave,1) - minvalue(wave,1);
            lowerthreshold = minvalue(wave, 1) + 0.2 .* height;
            higherthreshold = maxvalue(wave, 1) - 0.2 .* height;
            % Calculate the fall time
            for i = maxindex(wave, 1) : minindex(wave, 1)
                if (reversedata{channel}(wave, i) >= lowerthreshold) && (reversedata{channel}(wave, i + 1) < lowerthreshold)
                    endpoint = i + (lowerthreshold - reversedata{channel}(wave, i)) ./ (reversedata{channel}(wave, i + 1) - reversedata{channel}(wave, i));
                elseif (reversedata{channel}(wave, i) >= higherthreshold) && (reversedata{channel}(wave, i + 1) < higherthreshold)
                    startpoint = i + (higherthreshold - reversedata{channel}(wave, i)) ./ (reversedata{channel}(wave, i + 1) - reversedata{channel}(wave, i));
                end
            end
            falltime = endpoint - startpoint;
            % Calculate the previous rise time
            prevminnotfound = 1;
            prevminloc = maxindex(wave,1);
            while prevminnotfound
                if prevminloc > 1
                    if reversedata{channel}(wave, prevminloc) > reversedata{channel}(wave, prevminloc - 1)
                    
                        prevminloc = prevminloc - 1;
                    else
                        prevminnotfound = 0;
                    end
                else
                    prevminloc = 1;
                    prevminnotfound = 0;
                end
            end
            lowerthreshold = reversedata{channel}(wave, prevminloc) + 0.2 * (maxvalue(wave, 1) - reversedata{channel}(wave, prevminloc));
            higherthreshold = reversedata{channel}(wave, prevminloc) + 0.8 * (maxvalue(wave, 1) - reversedata{channel}(wave, prevminloc));
            for i = prevminloc : maxindex(wave, 1)
                if (reversedata{channel}(wave, i) <= lowerthreshold) && (reversedata{channel}(wave, i + 1) > lowerthreshold)
                    startpoint = i + (lowerthreshold - reversedata{channel}(wave, i)) ./ (reversedata{channel}(wave, i + 1) - reversedata{channel}(wave, i));
                elseif (reversedata{channel}(wave, i) <= higherthreshold) && (reversedata{channel}(wave, i + 1) > higherthreshold)
                    endpoint = i + (higherthreshold - reversedata{channel}(wave, i)) ./ (reversedata{channel}(wave, i + 1) - reversedata{channel}(wave, i));
                end
            end
            % Calculate the following rise time
            firstrisetime = endpoint - startpoint;
            follmaxnotfound = 1;
            follmaxloc = minindex(wave,1);
            while follmaxnotfound
                if follmaxloc < datapointnumber
                    if reversedata{channel}(wave, follmaxloc) < reversedata{channel}(wave, follmaxloc + 1)
                    
                        follmaxloc = follmaxloc + 1;
                    
                    else
                        follmaxnotfound = 0;
                    end
                else
                    follmaxloc = datapointnumber;
                    follmaxnotfound = 0;
                end
            end
            higherthreshold = minvalue(wave, 1) + 0.8 * (reversedata{channel}(wave, follmaxloc) - minvalue(wave, 1));
            lowerthreshold = minvalue(wave, 1) + 0.2 * (reversedata{channel}(wave, follmaxloc) - minvalue(wave, 1));
            for i = minindex(wave, 1) : follmaxloc
                if i == follmaxloc
                    endpoint = follmaxloc;
                    break
                end
                if (reversedata{channel}(wave, i) <= lowerthreshold) && (reversedata{channel}(wave, i + 1) > lowerthreshold)
                    startpoint = i + (lowerthreshold - reversedata{channel}(wave, i)) ./ (reversedata{channel}(wave, i + 1) - reversedata{channel}(wave, i));
                elseif (reversedata{channel}(wave, i) <= higherthreshold) && (reversedata{channel}(wave, i + 1) > higherthreshold)
                    endpoint = i + (higherthreshold - reversedata{channel}(wave, i)) ./ (reversedata{channel}(wave, i + 1) - reversedata{channel}(wave, i));
                end
            end
            secondrisetime = endpoint - startpoint;
            risefalltime{channel}(wave,1) = firstrisetime;
            risefalltime{channel}(wave,2) = 0 - falltime;
            risefalltime{channel}(wave,3) = secondrisetime;
        end
        %% calculate the width of the peak and trough
        troughinterval{channel}(wave,1) = channel;
        troughinterval{channel}(wave,2) = wave;
        peakinterval{channel}(wave,1) = channel;
        peakinterval{channel}(wave,2) = wave;
        threshold = troughthreshold(wave,1);
        troughstart = 0;
        troughend = 0;
        startfind = 0;
        endfind = 0;
        % Walking through the data to find the threshold crossing position
        % so as to find the trough position
        for i = 2: (datapointnumber - 1)
            if ((reversedata{channel}(wave,i - 1) > threshold)&&(reversedata{channel}(wave,i) <= threshold))
                troughstart = i - 1 + (threshold - reversedata{channel}(wave,i - 1)) / (reversedata{channel}(wave,i) - reversedata{channel}(wave,i - 1));
                startfind = 1;
            elseif((reversedata{channel}(wave,i - 1) <= threshold)&&(reversedata{channel}(wave,i) > threshold))
                troughend = i - 1 + (threshold - reversedata{channel}(wave,i - 1)) / (reversedata{channel}(wave,i) - reversedata{channel}(wave,i - 1));
                endfind = 1;
                if startfind
                    break
                end
            end
        end
        % Calculate the trough width
        if startfind && endfind
            troughinterval{channel}(wave,3) = troughstart;
            troughinterval{channel}(wave,4)= troughend;
            troughinterval{channel}(wave,5) = troughend - troughstart;  
        else
            troughinterval{channel}(wave,5) = 0;
        end
        % Calculate the trough area
        if troughinterval{channel}(wave,5)
            troughinterval{channel}(wave,6) = trapz(reversedata{channel}(wave, floor(troughstart):ceil(troughend)));
        else
            troughinterval{channel}(wave,6) = 0;
        end
        
        % Do the same process for a peak so as to find the peak width and
        % the peak area
        threshold = peakthreshold(wave,1);
        ppeakstart = 0;
        peakend = 0;
        startfind = 0;
        endfind = 0;
        for i = 2: (datapointnumber - 1)
            if ((reversedata{channel}(wave,i - 1) < threshold)&&(reversedata{channel}(wave,i) >= threshold))
                peakstart = i - 1 + (threshold - reversedata{channel}(wave,i - 1)) / (reversedata{channel}(wave,i) - reversedata{channel}(wave,i - 1));
                startfind = 1;
            elseif((reversedata{channel}(wave,i - 1) >= threshold)&&(reversedata{channel}(wave,i) < threshold))
                peakend = i - 1 + (threshold - reversedata{channel}(wave,i - 1)) / (reversedata{channel}(wave,i) - reversedata{channel}(wave,i - 1));
                endfind = 1;
                if startfind
                    break
                end
            end
        end
        if startfind && endfind
            peakinterval{channel}(wave, 3) = peakstart;
            peakinterval{channel}(wave, 4) = peakend;
            peakinterval{channel}(wave, 5) = peakend - peakstart;
        else
            peakinterval{channel}(wave,5) = 0;
        end
        if peakinterval{channel}(wave, 5)
            peakinterval{channel}(wave, 6) = trapz(reversedata{channel}(wave, floor(peakstart):ceil(peakend)));
        else
            peakinterval{channel}(wave, 6) = 0;
        end
            
        %% calculate fourier peak width
        FastFrourier{channel}(wave,:) = fft(reversedata{channel}(wave,:),NFFT)/L;
        FFABS{channel}(wave,:)=2*abs(FastFrourier{channel}(wave,1:NFFT/2+1));
        FFABSmean{channel}=mean(FFABS{channel},2);
        ffpeakinterval{channel}(wave,1) = channel;
        ffpeakinterval{channel}(wave,2) = wave;
        ffabsstd = std(FFABS{channel},0,2);
        threshold = FFABSmean{channel}(wave,1) + multiplier .* ffabsstd(wave,1);
        ffabsnumber = size(FFABS{channel},2);
        ppeakstart = 0;
        peakend = 0;
        startfind = 0;
        endfind = 0;
        % Walk through the FF transform result to find the threshold
        % crossing positions
        for i = 2: (ffabsnumber - 1)
            if ((FFABS{channel}(wave,i - 1) <= threshold)&&(FFABS{channel}(wave,i) > threshold))
                peakstart = i - 1 + (threshold - FFABS{channel}(wave,i - 1)) / (FFABS{channel}(wave,i) - FFABS{channel}(wave,i - 1));
                startfind = 1;
            elseif((FFABS{channel}(wave,i - 1) >= threshold)&&(FFABS{channel}(wave,i) < threshold))
                peakend = i - 1 + (threshold - FFABS{channel}(wave,i - 1)) / (FFABS{channel}(wave,i) - FFABS{channel}(wave,i - 1));
                endfind = 1;
                if startfind
                    break
                end
            end
        end
        if startfind && endfind
            ffpeakinterval{channel}(wave, 3) = peakstart;
            ffpeakinterval{channel}(wave, 4) = peakend;
            ffpeakinterval{channel}(wave, 5) = peakend - peakstart;
        else
            ffpeakinterval{channel}(wave, 5) = 0;
        end
        
        % column 1 of the property matrix: whether the midpoint of the wave
        % is the max point or not?
        if maxindex(wave,1) == sizewindow + 1
            waveformproperty{channel}(wave,1) = 1;
        else
            waveformproperty{channel}(wave,1) = 0;
        end
        % Column 2: position of the wave in the original file (used to
        % judge whether the classification result is right or not)
        waveformproperty{channel}(wave,2) = sorts.timestamp{channel}(1,wave);
    end
    [ffabsmax, ffabsmaxindex] = max(FFABS{channel},[],2);
    % 3: first rise/fall time related to the two maximum
    % 4: maximum interval's rise/fall time
    % 5: following rise/fall time
    % 6 trough width, 7 peak width, 8 power, 9
    % fast friorer power, 10 fast friorer peak width
    % 11 ff max value, 12 ff max position, 13 ff integral
    % 14 original wave integral, 15 max of the original wave, 16 min
    % of the original wave, 17 trough integral, 18 peak integral
    
    waveformproperty{channel}(:,3:5) = risefalltime{channel}(:,1:3);
    waveformproperty{channel}(:,6) = troughinterval{channel}(:,5);
    waveformproperty{channel}(:,7) = peakinterval{channel}(:,5);
    waveformproperty{channel}(:, 8) = standarddeviation;
    waveformproperty{channel}(:, 9) = ffabsstd;
    waveformproperty{channel}(:, 10) = ffpeakinterval{channel}(:, 5);
    waveformproperty{channel}(:, 11) = ffabsmax;
    waveformproperty{channel}(:, 12) = ffabsmaxindex(:, 1);
    waveformproperty{channel}(:,13) = trapz(FFABS{channel},2);
    waveformproperty{channel}(:,14) = trapz(reversedata{channel},2);
    waveformproperty{channel}(:,15) = maxvalue(:,1);
    waveformproperty{channel}(:,16) = minvalue(:,1);
    waveformproperty{channel}(:,17) = troughinterval{channel}(:,6);
    waveformproperty{channel}(:,18) = peakinterval{channel}(:,6);
end
timestamp = sorts.timestamp;
% noisermsmean will be used in later sections
noisermsmean = mean(rmsnoise);
clearvars -except file timestamp noisermsmean originalaverage FS T L NFFT f multiplier originaldata reversedata reversepositivedata originalpositivedata originalstd posoriginalstd sizewindow waveformproperty

%toc