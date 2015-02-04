% Generate data for plots
clc         % Clear work space
close all   % Close all figures
clear       % Clear all variables


% Input reaction parameters
reactionData = [1.5e3 3e-8 1.5e3 3e-8 100 5e-7 9.5e-5 1 1.5e-9 1.5e-9 2e-10	2e-9 1e-8 1e-8 500 5e6];

klpf = reactionData(1);
klmf = reactionData(2);
klpm = reactionData(3);
klmm = reactionData(4);
kap = reactionData(5);
kam = reactionData(6);
kn = reactionData(7);
cnmRatio = reactionData(8);
kbf = reactionData(9);
kbm = reactionData(10);  
y00 = reactionData(11);
p00 = reactionData(12);  
rs = reactionData(13);
rfm = reactionData(14); 
n = reactionData(15);
m = reactionData(16);


BasePath = 'C:\Users\youval\Documents\Education\Research\Prion line of defence\Paper 1\Plots\\Manuscrip plots\';
FilePath = [BasePath 'Figure8 - kinetic'];
FileName = 'PlotC.txt';



% Include parameters initialization code
eval('AggKineticsInitializationV16')
% Display data
eval('AggKineticsDisplayDataV16')

tic;

% Insert the main loop
eval('AggKineticsMainLoopV16')

rtime = toc;
disp(['Program run time: ' num2str(rtime)])
tMax = t(end);


disp (['treal [hours] = ' num2str(treal)])
disp ('+++++++++++++++++++++++++++++')
disp (['y(end) = ' num2str(y(end))])
disp (['g(end) = ' num2str(g(end))])
disp ('+++++++++++++++++++++++++++++')


%  Normelizing resualts - to avoid exponents in plots %
y1 = max(M);
y2 = max(F);
y3 = max(y0t);
y4 = max(g0t);

% ymax = max([y1 y2 y3 y4]);
ymax = max([y2 y3]);

yFactor = floor(log10(ymax));   % the 10 power of ymax
ymax = ymax/(10^yFactor);
M = M/(10^yFactor);
F = F/(10^yFactor);
y0t = y0t/(10^yFactor);
g0t = g0t/(10^yFactor);

% find good starting point to plot
% startPos is the time starting point
% eval('AggKineticsV16plotTimeStrat')
[tOnset,i] = AggKineticsV16plotTimeStrat(t);
% To cancel time sensative plotting set startPos = 1
startPos = 1;


% reducing number of points and changing the parameters names
t = t(startPos:10:end);            % Time
FABeta = y0t(startPos:10:end);    % Free A-Beta
F = F(startPos:10:end);            % Free aggregates
M = M(startPos:10:end);            % Membrane aggregates
g0 = g0t(startPos:10:end);         % Membrane PrPc monomers


figure()
h = gcf;
plot(t,FABeta,'--',t,F,'-',t,M,'-.');
legend('A-Beta','Free Agg.','Membrane Agg.','Location','NEO');
title('Total momomers in aggregates');
%plot label
xlabel(['Time  x10^{ ' num2str(xFactor) '} [Days]']);
ylabel(['Monomers concentration  x10^{ ' num2str(yFactor) '}[M]']);
% ymax = max([M F y0t]);
ylim([0 ymax*1.05]);
xlim([t(1) t(end)]);

% Save data to files
save([FilePath '\' FileName],'t','FABeta','F','M','g0','xFactor','yFactor');

