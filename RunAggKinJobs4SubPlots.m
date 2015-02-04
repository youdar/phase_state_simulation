% runing multipl jobs of AggKinetics
% Make sure AggKinetics is in fuction mode (That the function line is not
% commented out
% AggKineticsV9(klpf,klmf,klpm,klmm,kap,kam,kn,knm,kbf,kbm,y00,p00,rs,rfm,n,m,subPlotNum)

clear       % Clear all variables
clc         % Clear work space
close       % Close all figures


% Prodice set of 4 plots
% Vary yo to xo ratio

%                 (klpf,klmf,klpm,klmm,kap,kam ,kn  ,knm,kbf ,kbm  ,y00  ,p00,rs  ,rfm ,n  ,m  ,subPlotNum)
% AggKineticsV9Plots(3e4 ,2e-6,3e4 ,2e-6,1e2,2e-7,2e-5,2e7,2e-8,2e-10,1e-9,2e-9,1e-5,1e-8,500,1e7,1);
% clear       % Clear all variables
% clc         % Clear work space
% AggKineticsV9Plots(3e4 ,2e-6,3e4 ,2e-6,1e2,2e-7,2e-5,2e7,2e-8,2e-10,4e-10,2e-9,1e-5,1e-8,500,1e7,2);
% clear       % Clear all variables
% clc         % Clear work space
% AggKineticsV9Plots(3e4 ,2e-6,3e4 ,2e-6,1e2,2e-7,2e-5,2e7,2e-8,2e-10,1e-10,2e-9,1e-5,1e-8,500,1e7,4);
% clear       % Clear all variables
% clc         % Clear work space
% [xFactor,yFactor,xMax,yMax] = AggKineticsV9Plots(3e4 ,2e-6,3e4 ,2e-6,1e2,2e-7,2e-5,2e7,2e-8,2e-10,1.1e-10,2e-9,1e-5,1e-8,500,1e7,3);




% Vary a-beta flow - no sequestration
%                 (klpf,klmf,klpm,klmm,kap,kam ,kn  ,kbf ,kbm  ,y00  ,p00,rs,rfm ,n  ,m  ,subPlotNum)
AggKineticsV15Plots(3e4 ,2e-6,3e4 ,2e-6,1e2,2e-7,2e-5,2e-8,2e-10,2e-9,2e-9,0 ,3e-9,500,1e7,1);
clear       % Clear all variables
clc         % Clear work space
AggKineticsV15Plots(3e4 ,2e-6,3e4 ,2e-6,1e2,2e-7,2e-5,2e-8,2e-10,2e-9,2e-9,0 ,2e-8,500,1e7,2);
clear       % Clear all variables
clc         % Clear work space
AggKineticsV15Plots(3e4 ,2e-6,3e4 ,2e-6,1e2,2e-7,2e-5,2e-8,2e-10,2e-9,2e-9,0 ,5e-7,500,1e7,4);
clear       % Clear all variables
clc         % Clear work space
[xFactor,yFactor,xMax,yMax] = AggKineticsV15Plots(3e4 ,2e-6,3e4 ,2e-6,1e2,2e-7,2e-5,2e-8,2e-10,2e-9,2e-9,0 ,1e-7,500,1e7,3);





% Print label on plot
% xlabel(['Time  x10^{' num2str(xFactor) '} [H]'],'Position',[4.2,-1.65]);
% ylabel(['Number of monomers in aggregates  x10^{' num2str(yFactor) '}'],'Position',[-.5,0]);
xlabel(['Time  x10^{' num2str(xFactor) '} [H]'],'Position',[xMax*1.1,-0.12*yMax]);
ylabel(['Number of monomers in aggregates  x10^{' num2str(yFactor) '}'],'Position',[-0.12*xMax,yMax*1.1]);
