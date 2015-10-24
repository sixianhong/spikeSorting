% Find spikes from the collected wave data
% Most of the code come from Dr. Langhals's code
% Modified line 372 to 395 to denoise the signal using Reverse PCA
% Store the denoised signal in reversedata, and the original signal in originaldata
% Deleted the clustering part from Dr. Langhals's code

clear all; close all; clc;

% Start a pool of workers if one is not already running
%{
 if ~matlabpool('size')
     matlabpool open
 end
%}

% Default options for variable name of spike data.  If these are in the
% file, it doesn't prompt you and just uses them to analyze.
spikedatafields = {'STRM', 'fSig', 'DDT', 'PLX', 'TDT'};

% Undo RTCAR & add offline CAR (1 = undo, 0 = leave)
CARundo = 1;
CARadd = 0;

% Noise multiplier level (i.e. 6*SD for noise)
stdmult = 6;

% Theshold for spike detection in SD's
detthresh = 3.5;
snrthresh = 1;
fcmthresh = 0.8;

% Option to detect spikes with a large positive deflection with
% insufficient negative to cross threshold.  Set to 1 to detect.
positivespikes = 1;

%The current program is set so that maxunits must be at least 3.
maxunits = 7;

%Set this to 1 to plot all clusters, or 0 to only plot counted ones.
allclusters = 0;

% Set default filter cutoffs
lfplowfreq = 1;
lfphighfreq = 100;
spikelowfreq = 300;
spikehighfreq = 5000;

% Amount to decimate LFP data from original sampling rate for plotting
lfpdec = 10;

% Spike window duration in seconds
window = .003;

% Time to start using the data (eliminate transients from turnon)
starttime = 0;

% Initialize Figure Properties (set savefigs to 1 if save, 0 if not).
savefig = 0;
savetif = 1;
saveeps = 0;
set(0,'DefaultFigureWindowStyle','docked');
rect = [0.25, 0.25, 10.5, 8];
plotColor = ['b','g','r','c','m','y','k'];
plotMarker = ['o','x','s','d','*','p','h','+'];

% Initialize some variables
spikedatafield = [];
intscale = 1e7;

% Pick the file to analyze
[fnames, file.pathname] = uigetfile('*.mat', 'Select Files', 'MultiSelect', 'on');
if file.pathname
    cd(file.pathname)
    if iscell(fnames)
        file.fname = sort(fnames);
    else
        file.fname = sort({fnames});
    end
    
else
    disp('No File Selected')
    break
end
clear fnames

% Check for Analysis directory or create one
if ~isdir('Analysis')
    mkdir('Analysis')
end
file.savepath = [file.pathname 'Analysis'];

tic
for multifile = 1:numel(file.fname)

    % Clear previous run figures and data if present
    clear data %Here data is the variable name, maybe. All the data files saved are named data
    close all
    
    load([file.pathname file.fname{multifile}])
    dfields = fieldnames(data); %All the variable names in the data file
    if ~isfield(data, spikedatafield)
        % spikedatafields contains preset options for streamed spike data
        spikedataindexes = isfield(data, spikedatafields);
        
        % Determine field containing spike data with error checking
        if sum(spikedataindexes) == 1
            spikedatafield = char(spikedatafields(spikedataindexes));
        elseif numel(dfields) == 1
            spikedatafield = dfields{1};
        else
            if sum(spikedataindexes)
                temp = spikedatafields(spikedataindexes);
            elseif isfield(data, 'eventnames')
                temp = data.eventnames;
            else
                temp = dfields;
            end
            [tempselect,ok] = listdlg('PromptString','Select the Variable with you data','ListString',temp, 'SelectionMode', 'single');
            if ok
                spikedatafield = temp{tempselect};
            else
                h = warndlg('No Variable selected');
                close(1);
                return
            end
            clear ok tempselect temp
        end
    end
    
    % Convert integer data to float if necessary on the first iteration of
    % the loop.  Does not support mixed int and float in same analysis
    % bunch.
    if multifile == 1
        if isinteger(data.(spikedatafield).data{1})
            scalebutton = questdlg('Is the scale factor of your integer data 1e7?','Integer Conversion Factor');
            switch scalebutton
                case 'Yes'
                    intscale = 1e7;
                case 'No'
                    tempans = inputdlg('Enter your scale factor - Must be >1');
                    intscale = str2double(tempans);
                    if isnan(intscale)
                        disp('Invalid Number - Canceling Program')
                        break
                    end
                case 'Cancel'
                    disp('Canceling Program')
                    break
            end
        else
            % Sets this to 1 to not affect floating data
            intscale = 1;
        end
    end
    
    % Determine number of channels in data
    numchan = numel(data.(spikedatafield).data);
    
    %Round sampling rate
    %     samprate = round(data.(spikedatafield).samprate(1));
    %samprate = round(data.(spikedatafield).samprate);
    samprate=25000;
    
    % Change endsample to use limited amount of data if too large.
    %     endsample = min(cellfun(@numel, data.(spikedatafield).data));
    endsample = cellfun(@numel, data.(spikedatafield).data);                        %The total amount of data points in each channel
    
    % Add RTCAR signal back into original
    if isfield(data, 'CARs')
        if CARundo == 1
            % minsize is used to make sure both CAR and STRM are same size
            minsize = min(length(data.CARs.data{1}), endsample);
            % NOTE: THIS NEEDS TO BE FIXED! Resaving onto same structure
            % can temporarily multiply the size by number of channels using
            % this cellfun method.
            data.(spikedatafield).data = cellfun(@(x) x(1:minsize)+data.CARs.data{1}(1:minsize), data.(spikedatafield).data, 'UniformOutput', false);
            if CARadd == 1
                CARname = 'ppCAR';
            else
                CARname = 'noCAR';
            end
        elseif CARadd == 1
            disp('You cant add CAR twice');
            break;
        else
            CARname = 'rtCAR';
        end
    else
        CARname = '';
    end
    
    

    % Do CAR Processing
    if CARadd == 1
        cardata = zeros(1,min(endsample));
        CARname = 'CAR';
        for chan = 1:numchan
            cardata = cardata+data.(spikedatafield).data{chan}/numchan;
        end
    end
    
    % THIS IS THE POINT WHERE WE ARE GOING TO DO SINGLE CHANNEL PROCESSING!
    % All initializations moved to above the channel for loop
    
    % grabs and filters data
    datafilter = cell(1,numchan);
    datafilterlfp = cell(1,numchan);
    snip = cell(numchan,1);
    chanminholder = cell(numchan,1);
    sorts.timestamp = cell(numchan,1);
    rmsnoise = zeros(1,numchan);
    ppnoiseraw = zeros(1,numchan);
    ppnoiselfp = zeros(1,numchan);
    sorts.snr = zeros(numchan,1);
    correlationmatrix = [];
    correlationmatrix2 = [];
    sorts.meansig = zeros(numchan,1);
    ppnoise = zeros(numchan,1);
    sorts.meanwaveforms = cell(numchan,1);
    sorts.snip = cell(numchan,1);
    sorts.sniptime = cell(numchan,1);
    U = cell(numchan,1);
    maxU = cell(numchan,1);
    fcmdata = cell(numchan,1);
    reversedata=[];
    originaldata=[];
    
    
    
    squareplot = ceil(sqrt(numchan));
    
    % hh = waitbar(0,'Analyzing Channels...');
    
    for channel = 1:numchan
        % If there's only 1 samplerate per channel, use the first one on
        % all channels
        clear X;
        if channel>numel(samprate)
            samprate(channel) = samprate(1);
        end
        
        %Determine starting sample point after delay
        startsample = starttime*samprate(channel);
        
        % Continue to next for loop iteration if there is no data in that
        % channel or only LFP data
        if ~endsample(channel)||samprate(channel)<10000||startsample>endsample(channel)
            continue
        end
        
        sampratelfp = round(samprate(channel)/lfpdec);

        % Calculate filter coefficients
        lfpfilter = [max(0.001,lfplowfreq/samprate(channel)*0.5) min(0.999,lfphighfreq/(samprate(channel)*0.5))];
        [blfp, alfp] = butter(1, lfpfilter);
        % If using lfp data, you need to limit the filter so that the
        % butterwoth doesn't have invalid filter points based on the sample
        % rate
        spikefilter = [max(0.001,spikelowfreq/samprate(channel)*0.5) min(0.999,spikehighfreq/(samprate(channel)*0.5))];
        [b, a] = butter(1, spikefilter);
        
        % Determine timing for data
        sizewindow = ceil(window/2*samprate(channel));    %window in data points
        spiketimeaxis = 1000*(0:1:(2*sizewindow-1))/samprate(channel); %
        
        seconds = numel(data.(spikedatafield).data{1})/samprate(channel)-starttime;
        if seconds<0
            disp('Less than "starttime" of data in file - unable to analyze')
            continue
        end


        datafilter = filter(b,a,double(data.(spikedatafield).data{channel}(:, startsample+1:endsample(channel)))/intscale, [],2);
        % VERIFY DECIMATE is not overly affecting data - it incorporates an
        % 8th order lowpass before resampling.  NOTE: this is on top of the
        % bandpass done below using a butterworth filter.
        datafilterlfp = decimate(filter(blfp,alfp,double(data.(spikedatafield).data{channel}(:, startsample+1:endsample(channel)))/intscale, [],2), lfpdec);
        
        % calculates the rms of each of the 16 channels
        rmsnoise(channel) = norm(datafilter)/sqrt(numel(datafilter));
        
        ppnoiseraw(channel) =stdmult*std(datafilter); %determine the cut-off voltage
        ppnoiselfp(channel) = stdmult*std(datafilterlfp);
        originalaverage{channel} = mean(datafilter);
        originalstd{channel} = std(datafilter);
        
        pthresh = mean(datafilter) + detthresh*std(datafilter);
        nthresh = mean(datafilter) - detthresh*std(datafilter);
        
        holdermax = find(datafilter > pthresh);
        holdermin = find(datafilter < nthresh);
        
        chantemp = [];
        chantempholder = [];
        
        placeholder = 0;
        
               %sample waveform is negative peak first
        % Go through all crossings
        for i = 1:numel(holdermin)
            % Is it too early?
            if holdermin(i) <= (sizewindow)+1
                % And is it too close to the end?
            elseif holdermin(i) <= (numel(datafilter) - sizewindow)
                % Figure out minimum value and point at which it occurs in
                % snippet window
                [mintemp, mintempposition] = min(datafilter(holdermin(i)-sizewindow:holdermin(i)+sizewindow-1));
                % Adjust minimum position to be the index for the whole
                % array, not just the small snippet fed into it.
                mintempposition = holdermin(i) - sizewindow-1 + mintempposition;
                % If this minimum happens to be the holder value in all
                % minimums, then continue
                if (mintemp == datafilter(holdermin(i)));
                    placeholder = placeholder +1;
                    chantemp(placeholder, :) =  datafilter(mintempposition-sizewindow:mintempposition+sizewindow-1);
                    chantempholder(size(chantempholder,2)+1:size(chantempholder,2)+2*sizewindow) = (mintempposition-sizewindow):(mintempposition+sizewindow-1);
                    sorts.timestamp{channel}(placeholder) = (mintempposition+startsample);
                end
            end
        end
        
        if positivespikes == 1
            % sample waveform is positive peak first
            for i = 1:length(holdermax)
                if holdermax(i) <= (sizewindow)+1
                elseif holdermax(i) <= (numel(datafilter) - sizewindow)
                    [maxtemp, maxtempposition] = max(datafilter(holdermax(i)-sizewindow:holdermax(i)+sizewindow-1));
                    [mintemp, mintempposition] = min(datafilter(holdermax(i)-sizewindow:holdermax(i)+sizewindow-1));
                    maxtempposition = holdermax(i) - sizewindow-1 + maxtempposition;
                    if (maxtemp == datafilter(holdermax(i))) && maxtemp>=abs(mintemp) && mintemp>nthresh
                        placeholder = placeholder +1;
                        chantemp(placeholder, :) =  datafilter(maxtempposition-(sizewindow):maxtempposition+(sizewindow-1));
                        chantempholder(size(chantempholder,2)+1:size(chantempholder,2)+2*sizewindow) = (maxtempposition-(sizewindow)):(maxtempposition+(sizewindow-1));
                        sorts.timestamp{channel}(placeholder) = (maxtempposition+startsample);
                    end
                end
            end
        end
        
        % Verify that candidates were found, if so, save in snip
        if size(chantemp,1)>0
            snip{channel}(1:placeholder, :) = chantemp;
            chanminholder{channel}(1:length(chantempholder)) = chantempholder;
        else
            % Setup conditions for set diff below for situations with no
            % spikes
            snip{channel}(1, 1:2*sizewindow+1) = zeros(1,2*sizewindow+1);
            chanminholder{channel}(1) = 0;
        end
        
        clear chantemp
        
        %subtracts out signal and calculates the noise for each channel
        indextemp = 1:numel(datafilter);
        % Determine if spikes were in holder
        nosignal1 = find(chanminholder{channel} > 0);
        % Grab non-zero indexes
        nosignal2 = chanminholder{channel}(nosignal1);
        % Compare indexes with spike data and grab non-spike indexes
        nosignal3 = setdiff(indextemp,nosignal2);
        % Grab all non-spike data
        nosignal4 = datafilter(nosignal3);
        % Calculate noise on the channel
        ppnoise(channel) =  stdmult*std(nosignal4);
        
        %holder checks to see there is a significant number of points above
        %on the given channel above threshold;  If there isn't enough points to
        %make reasonable clusters, the cluster is dropped from further analysis
        holder = find(sum(snip{channel},2));
        if length(holder) > 30
            
            %5 dimensional PCA
            [PC, SCORE, LATENT, TSQUARE] = pca(squeeze(snip{channel}(1:holder(length(holder)),:)));
            X = SCORE;
            columnofX=size(X,2);
            rowofX=size(X,1);
            for i=1:rowofX
                A(i,1)=0;
            end
            for i=5:columnofX
                X(:,i)=A;
            end
            clear A;
            newdata=(PC)*(X');
            olddata=(PC)*(SCORE');
            initialdatamean=mean(snip{channel}(1:holder(length(holder)),:),1);
            newdata=newdata';
            for i=1:rowofX
                newdata(i,:)=newdata(i,:)+initialdatamean;
            end
            olddata=olddata';
            for i=1:rowofX
                olddata(i,:)=olddata(i,:)+initialdatamean;
            end
            reversedata{channel}=newdata;
            originaldata{channel}=olddata;
            
        end
        
        
    end
    
end

clearvars -except file sorts originalstd rmsnoise samprate originalaverage originaldata reversedata reversepositivedata originalpositivedata originalstd posoriginalstd sizewindow