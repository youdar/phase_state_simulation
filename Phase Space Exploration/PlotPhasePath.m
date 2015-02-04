clc         % Clear work space
close all   % Close all figures
clear       % Clear all variables
% Plot phase path
FileName = 'PhaseLine vary y0 and rfm-rs.txt';
% FileName = 'PhaseLine vary y0 and p0.txt';
% FilePath = 'C:\Users\youval\Documents\Education\Research\Prion line of defence\Paper 1\Matlab workspace\Critical aggregation phase line\';
% FilePath = 'C:\Users\youval\Documents\Education\Research\Prion line of defence\Paper 1\Matlab workspace\40% a-beta reduction phase line\';
FilePath = 'C:\Users\youval\Documents\Education\Research\Prion line of defence\Paper 1\Matlab workspace\40% a-beta reduction phase line\Test runs\';
% select number of base parameter
basePnum = 11;
% select adjustable parameter
adjPnum = 14;
SemiLogPlot = 1;    % 0: regular plot, 1: semilog plot
% Choose which parameters to fix:
% 1: klpf       elongation (monomer attachment) for free aggregates
% 2: klmf       monomer dettachment, for free aggregates 
% 3: klpm       elongation (monomer attachment) for membrane bounnd aggregates
% 4: klmm       monomer dettachment, for membrane bounnd aggregates
% 5: kap        attachemnt rate, of aggregates to te membrane
% 6: kam        dettachment rate, of aggregates to te membrane
% 7: kn         free aggregates nucleation rate
% 8: cnmRatio   cnmRatio = (cnm equivalent)/cn = cnm*p0/cn
% 9: kbf        free aggregates brakage rate
% 10: kbm       mambrane aggregates brakage rate
% 11: y00       initial A-Beta concentration
% 12: p00       initial Prion concentration
% 13: rs        Sequestration rate
% 14: rfm       Out-flow rate
% 15: n         size of maximum aggregates length
% 16: m         number of iterations

% nFilePath = input(['Current path is ' FilePath '  : '],'s');
% if ~isempty(nFilePath)
%     FilePath = nFilePath;
% end
% nFileName = input(['Current file name is ' FileName '  : '],'s');
% if ~isempty(nFileName)
%     FileName = nFileName;
% end

fileID = fopen([FilePath FileName],'r');       % Open file for reading
if fileID<0
    disp('_______________________________________________________________')
    disp('Error !!! Can not open data file !!!')
    disp('Please check file name and path')
    disp([FilePath FileName])
    disp('_______________________________________________________________')
end

ParameterArray = zeros(1,16);    % Set the parameter array size
ParametersName = char('{\itk}_{lf}^+','{\itk}_{lf}^-','{\itk}_{lm}^+','{\itk}_{lm}^-','{\itk}_a^+','kam','kn','cnmRatio','{\itk}_{bf}','{\itk}_{bf}','A\beta','PrPc','rs','rfm');



y = zeros(1000,1);              % Initializing point arrays
x = zeros(1000,1);              % Initializing point arrays


DataLine = fgetl(fileID);
i = 1;
while ischar(DataLine)
    if DataLine(1)=='['          % Collecting data only from lines that starts with [
        ParameterArray = eval(DataLine);
        y(i) = ParameterArray(basePnum);      % Set initial value
        x(i) = ParameterArray(adjPnum);       % Set initial value   
        i = i+1;
    end
%     disp(mat2str(DataLine))
    DataLine = fgetl(fileID);
end

fclose(fileID);
Endi = max(find(x,1,'last'),find(y,1,'last'));
% taking only data points
x = x(1:Endi);
y = y(1:Endi);
% sorting the points
xy  = sortrows([x y]);
x = xy(:,1);
y = xy(:,2);

% x = x(1:end-1);
% y = y(1:end-1);

% midCol = round(Endi/2);
% midCol = Endi - midCol;
% % x = x(midCol:end);
% % y = y(midCol:end);
% x = x(1:midCol);
% y = y(1:midCol);

% [xFactor,yFactor] = XYfactor(x(end),y(end));
[xFactor,yFactor] = XYfactor(x(1),y(1));

x = x/xFactor;
y = y/yFactor;
[fittedFun,fitParameter] = FitFunctionDataPlot(x,y);
% fittedFun
% fitParameter
disp(fittedFun)
disp('Fit Parameters'); disp(fitParameter)
FitFunctions = ['    poly1     ','poly2     ','poly3     ','weibull    ','fourier4    ','rat11    ','rat22     ','power2     ','exp1     ','exp2'];
disp(FitFunctions)
% plot label
% xlabel([ParametersName(adjPnum) ' x10^{ ' num2str(floor(log10(xFactor))) '}']);
% ylabel([ParametersName(basePnum) ' x10^{ ' num2str(floor(log10(yFactor))) '}']);
xpowerstr = num2str(floor(log10(xFactor)));
ypowerstr = num2str(floor(log10(yFactor)));
xnamestr = ParametersName(adjPnum,:);
ynamestr = ParametersName(basePnum,:);

h = figure();
xx = x*xFactor;
y2 = feval(fittedFun,x);
if SemiLogPlot
    semilogx(xx,y,'o',xx,y2,'r')
    hold on
    textlength = find(isletter(xnamestr),1,'last');
    xnamestr = xnamestr(1:textlength);
    strx = ['Log[' xnamestr ']'];
else
    plot(xx,y,'o',xx,y2,'r')
    hold on
    strx = [xnamestr ' [\times 10^{ ' xpowerstr '}]'];
end



stry = [ynamestr '[\times 10^{ ' ypowerstr '}] M'];
xlabel(strx,'FontSize',12);
ylabel(stry,'FontSize',12);
% title(FileName(1:end-4))
% if abs(y-1)<0.1
%     ylim([0.5 1.5])
% end
% ylim([1.5 1.8])
% Saving the plot to a file
filename = [FilePath FileName(1:end-4) ' data fit  ' '.jpg'];
% saveas(h,filename);
hold off