% Makes Fig. S1: example for time of observation bias

time1 = 1:20;
time2 = 20:48;
clf 
T1 = 3*sin(2*pi*(time1 - 9)/24);
T2 = 1.5*sin(2*pi*(time2- 9)/24);
d = T2(1) - T1(end);
T2 = T2 - d;

T = smoothPH([T1 T2(2:end)],3);
time = 1:48;

timePlot = 1:1/60:48;
Tplot = interp1(time,T,timePlot,'spline');

plot(timePlot, Tplot, 'k','linewidth',2)

hold on 
plot([7 7],[-3 3],':k')
plot([17 17],[-3 3] ,'--k')
plot([31 31],[-3 3],':k')
plot([41 41],[-3 3],'--k')

xlim([1 48])
set(gca,'fontsize',10)
set(gca,'box','on')
set(gca,'xtick',[4:4:48])
set(gca,'xticklabel',[4:4:24 4:4:24])

set(gcf, 'Color', 'w');
figName = [figDir '/timeofObsExample.pdf'];
export_fig(figName, '-m5','-a1')