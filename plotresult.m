% Plot step
channels = size(newcatgory,2);
ro = size(newcatgory{1}{1},1);
a = 10;

for chn = 1:a
    h(chn) = figure(chn);
    set(gcf, 'name', sprintf('Channel %d', chn));
    groupnum = size(newcatgory{chn},2);
    for group = 1: groupnum
        figure(chn)
        wavenum = size(newcatgory{chn}{group},1);
        subplot(2,3,group)
        datapointnum = size(newcatgory{chn}{group},2);
        ro = 1:datapointnum;
        for i = 1: wavenum
            %plot(ro, newcatgory{chn}{group}(i,:));
            plot(ro, neworiginalshape{chn}{group}(i,:));
            hold on;
            title(sprintf('Cluster %d', group));
            xlabel('Time point')
            ylabel('Amplitude(Normalized)')
            axis([0 76 -1.5 1.5]);
            %axis([0 45 -7e-5 7e-5])
        end
    end
    cd(file.savepath)
    print(h(chn), '-djpeg', [get(h(chn), 'name')])
end
% Plot out the original shape
figure
r = 1:size(comparison.shape,2);
for i = 1:10
    subplot(3,5,i)
    plot(r,comparison.shape(i,:));
    axis([0 38 -1.5 1.5])
end