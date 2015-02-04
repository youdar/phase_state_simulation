function [yLast,gLast,tLast,TestFlag,TestFlagValue,n10,n90,x1,x2,alarm,tt,ft,F,M]=AggKineticsV16PhaseExplor(DataArray)
%
% ft is the quantity that is being tested. change it at the end of the
% function
%

    
% change random numbers stream
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));

% get data from the DataArray

klpf = DataArray(1); 
klmf = DataArray(2);
klpm = DataArray(3);
klmm = DataArray(4);
kap = DataArray(5);
kam = DataArray(6);
kn = DataArray(7);
cnmRatio = DataArray(8);        % cnmRatio = (cnm equivalent)/cn = cnm*p0/cn
kbf = DataArray(9);
kbm = DataArray(10);
y00 = DataArray(11);
p00 = DataArray(12);
rs = DataArray(13);
rfm = DataArray(14);
n = DataArray(15);
m = DataArray(16);

% Normalization factor
% Don't forget to adjust real time accordingly
alarm = 0;

% Include parameters initialization code
eval('AggKineticsInitializationV16')
% Insert the main loop
eval('AggKineticsMainLoopV16')


if y00<0 || p00<0
    TestFlag = 1;                   % Test fails -> TestFlag=1
    TestFlagValue = 1;
else
    % Test how many of the last 500 steps had segnificant aggregation    
    %  TestVal = 100*F(end-500:end)/(F(end-500:end)+y0t(end-500:end));  
    % for reduced A-Beta use 
    [TestFlag,TestFlagValue,n10,n90,x1,x2] = TestABetaValue(y0t,t);
    % For critical aggregation use
%     [TestFlag,TestFlagValue,n10,n90,x1,x2] = TestAggValue(F,M,t);
    %  I take out the last term since in high flow rates it can cause
    %  problems 
end

 
% ft = F;
ft = y0t;
% convert concentration
F = F*vScale/v;
M = M*vScale/v;
% y0t = y0t*vScale/v;
% g0t = g0t*vScale/v;

% Reduce array size
tt = [t(1:9) t(10:10:end)]*(10^xFactor);    % Time in days
ft = [ft(1:9) ft(10:10:end)];
n10 = round((n10 - 9)/10);
n90 = round((n90 - 9)/10);
% n90 = find(tt>t(n90),1);
% n10 = find(tt>t(n10),1);

yLast = y(end);
gLast = g(end);
tLast = tt(end);
end