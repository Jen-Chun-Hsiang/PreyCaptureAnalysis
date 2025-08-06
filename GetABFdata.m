addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BasicFunctions\OpenSource');
filename = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\ephys\c120723_0003.abf';
d = abfload(filename,'start',1,'stop','e');
%%
d = data;
figure; 
for i = 1:3
    subplot(3, 1, i);
    maxv = max(d(:, i));
    minv = min(d(:, i));
    plot(d(:, i));
    ylim([minv maxv*1.2])
    box off
end
