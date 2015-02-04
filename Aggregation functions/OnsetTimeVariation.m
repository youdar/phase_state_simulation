% collect information on onset time variation

clc         % Clear work space
close all   % Close all figures
clear       % Clear all variables
    
    
N = 50;
onsettime = zeros(1,N);
for ii=1:N
    % Input reaction parameters
    eval('AggKineticsV16Input')
    % Include parameters initialization code
    eval('AggKineticsInitializationV16')
    % Display data
%     eval('AggKineticsDisplayDataV16')
%     disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
%     disp(['Iteration #: ' num2str(ii)])
    % Insert the main loop
    eval('AggKineticsMainLoopV16')
    %  Normelizing resualts - to avoid exponents in plots %
%     y1 = max(M);
%     y2 = max(F);
%     y3 = max(y0t);
%     y4 = max(g0t);
%     ymax = max([y2 y3]);
%     yFactor = floor(log10(ymax));   % the 10 power of ymax
%     ymax = ymax/(10^yFactor);
%     M = M/(10^yFactor);
%     F = F/(10^yFactor);
%     y0t = y0t/(10^yFactor);
%     g0t = g0t/(10^yFactor);
    [tj,j] = AggKineticsV16plotTimeStrat(t);
%     figure(1)
%     plot(t,y0t,'--',t,F,'-',t(1:10:end),M(1:10:end),'-.',[tj,tj],[0,ymax]);
%     ttlstr = ['Iteration #: ' num2str(ii)];
%     title(ttlstr)
%     drawnow
    onsettime(ii) = tj*(10^xFactor);
end




