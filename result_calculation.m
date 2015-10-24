% Read in the comparison matrix that is used to judge the accuracy of the
% algorithm

load([file.pathname 'comparison.mat']);

% Plot out the original shape
figure
r = 1:size(comparison.shape,2);
for i = 1:10
    subplot(3,5,i)
    plot(r,comparison.shape(i,:));
    axis([0 38 -1.5 1.5])
end
clear chn a group classificationresult datapointnum ro wavenum targetgroup targettimepoint pointsfound range ifininterval num position pointsfound foundposition
% Calculate the accuracy result
channel = size(newproperty,2);
datapointnum = size(newcatgory{1}{1},2);
for chn = 1:channel
    totalgroup = size(newcatgory{chn},2);
    shapetemplate = find(comparison.template(chn,:));
    for group = 1:totalgroup
        wavenum = size(newcatgory{chn}{group},1);
        for templatenum = 1:size(shapetemplate,2)
            % Find the time position for wave to occur
            indicies = find(comparison.time(shapetemplate(1,templatenum),:)>0);
            targettimepoint = comparison.time(shapetemplate(1,templatenum),indicies);
            pointsfound = 0;
            for i = 1: wavenum
                range(i,1) = newproperty{chn}{group}(i,2) - datapointnum/2;
                range(i,2) = newproperty{chn}{group}(i,2) + datapointnum/2;
                ifininterval = (targettimepoint>=range(i,1)) & (targettimepoint<=range(i,2));
                [num,position] = find((targettimepoint>=range(i,1)) & (targettimepoint<=range(i,2)) > 0);
                if ~isempty(position)
                    pointsfound = pointsfound + 1;
                    foundposition{pointsfound} = position;
                end
            end
            expectedwavenum = size(targettimepoint,2);
            foundwavenum = pointsfound;
            classifiedwave = wavenum;
            % row means which template, column 1 is which is the template
            % column 2 compares the sorting
            % result with the expected "correct answer", column 3 compares the
            % amoung correctly sorted to the total amount being classified
            % into the group
            classificationresult{chn}{group}(templatenum, 1) = shapetemplate(templatenum);
            classificationresult{chn}{group}(templatenum, 2) = foundwavenum / expectedwavenum;
            classificationresult{chn}{group}(templatenum, 3) = foundwavenum / classifiedwave;
        end
    end
end
for chn = 1:channel
    totalgroup = size(newcatgory{chn},2);
    for group = 1:totalgroup
        sortedresult{chn}{group} = sortrows(classificationresult{chn}{group},2);
        totaltemplate = size(sortedresult{chn}{group},1);
        bestclassificationresult{chn}(group,:) = sortedresult{chn}{group}(totaltemplate,:);
    end
end