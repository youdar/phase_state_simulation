% Prodice set of 4 plots


clear       % Clear all variables
clc         % Clear work space
close       % Close all figures
% change random numbers stream
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));

% Start labs
matlabpool                      % allocate local MATLAB labs for parallel functions execution  
labSize = matlabpool('size');   % find out how many processors we use
if labSize~=4                   % test if we have enough labs
    matlabpool close
    disp(['Only ' num2str(labSize) ' are available !!!!'])
    disp('This program is designed for 4 labs')
    JobStatus = 0;              % Indicate that we have Error
else
    disp(['Doing ' num2str(labSize) ' runs at paralell'])
    disp('          ')
end

% update and prepare the inputMatrix matrix for the 4 labs
%                  [klpf,klmf,klpm,klmm,kap,kam     ,kn  ,kbf  ,kbm  ,y00      ,p00 ,rs       ,rfm      ,n  ,m]
inputMatrix(1,:) = [3e4 ,2e-6,3e4 ,2e-6,1e2,2.042e-7,2e-5,2e-10,2e-10,1.68e-010,2e-9,2.8e-010 ,1.16e-011,500,1.2e7];
inputMatrix(2,:) = [3e4 ,2e-6,3e4 ,2e-6,1e2,2.042e-7,2e-5,2e-10,2e-10,1.68e-010,2e-9,2.8e-010 ,1.16e-010,500,2.1e7];    % rfm
inputMatrix(3,:) = [3e4 ,2e-6,3e4 ,2e-6,1e2,2.042e-7,2e-5,2e-10,2e-10,1.68e-010,  0 ,2.8e-010 ,1.16e-011,500,1.5e7];    % p00
inputMatrix(4,:) = [3e4 ,2e-6,3e4 ,2e-6,1e2,2.042e-7,2e-5,2e-10,2e-10,1.68e-010,  0 ,2.8e-010 ,1.16e-010,500,2.3e7];    % p00,rfm

tic;
spmd
    [xFactor,yFactor,xMax,yMax,t,M,F,y0t,g0t] = AggKineticsV10Plots(inputMatrix(labindex,:));
end
toc

h = figure();
j = [1 2 4 3];
for ii=1:4           % i is the subPlotNum
    i = j(ii);
    switch i
        case 1
            subplot('position',[.09 .56 .4 .4])
        case 2
            subplot('position',[.56 .56 .4 .4])
        case 3
            subplot('position',[.09 .09 .4 .4])
        case 4
            subplot('position',[.56 .09 .4 .4])
    end
    tt = t{i};
    plot(tt,M{i},'--',tt,F{i},'-',tt,y0t{i},'-.',tt,g0t{i},':');
    ytop = yMax{i}*1.05;
    xtop = tt(end);
    ylim([0 ytop]);
    xlim([0 xtop]);

    % do this after the plotting
    switch(i)
        case 1
%             set(gca,'XTick',[])
%             ylabel('[#]');
            xlabel(['x10^{' num2str(xFactor{i}) '}'],'Position',[xMax{i}*1.02,0*yMax{i}],'FontSize',8);
            text(xtop*0.94,0.94*ytop,'A');
        case 2
%             set(gca,'XTick',[])
            set(gca,'YTick',[])
            xlabel(['x10^{' num2str(xFactor{i}) '}'],'Position',[xMax{i}*1.02,0*yMax{i}],'FontSize',8);
            text(xtop*0.94,0.94*ytop,'B');
        case 3
%             xlabel('t[h]');
%             xlabel('t');
%             ylabel('[#]');
            xlabel(['x10^{' num2str(xFactor{i}) '}'],'Position',[xMax{i}*1.02,0*yMax{i}],'FontSize',8);
            text(xtop*0.94,0.94*ytop,'C');
        case 4
            set(gca,'YTick',[])
            xlabel(['x10^{' num2str(xFactor{i}) '}'],'Position',[xMax{i}*1.02,0*yMax{i}],'FontSize',8);
            text(xtop*0.94,0.94*ytop,'D');
%             xlabel('t[h]');
%             xlabel('t');
    end
end


% Print label on plot
% xlabel(['Time  x10^{' num2str(xFactor) '} [H]'],'Position',[4.2,-1.65]);
% ylabel(['Number of monomers in aggregates  x10^{' num2str(yFactor) '}'],'Position',[-.5,0]);
i = 3;
text(xMax{i}*0.95,-0.16*yMax{i},'Time in [H]');
ylabel(['Number of monomers in aggregates  x10^{' num2str(yFactor{i}) '}'],'Position',[-0.12*xMax{i},yMax{i}*1.1]);

filename = ['Four plots - AB aggregation' '.fig'];
saveas(h,filename);
