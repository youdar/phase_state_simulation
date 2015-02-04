% This program is used to explore a critical trajectory for two paramters
% From our set of 16 parameters:
% {klpf,klmf,klpm,klmm,kap,kam,kn,cnmRatio,kbf,kbm,y00,p00,rs,rfm,n,m,x}
% we are fixing all but two.
% then choose one of the two as the base paramerter and the other as the adjustable one.
% We choose a base value and an adjustable value on the critical or desire phase space line
% then vary the base value and adjust the other variable until we return
% to the desire phase line.

% Clean work space
clc             % Clear work space
close all       % Close figures
clear           % Clear all variables
% change random numbers stream
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));

% Initial time scale settings
% tDay = 60*60*24;    % Seconds in a day
% tHour = 3600;       % Seconds in one hour
% % tFactor = 1e-6;     % Time factor for simulation
% tFactor = 1;
% rs = tFactor/tHour;
% rfm = tFactor/tDay;

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
% 14: rfm       Clearance rate
% 15: n         size of maximum aggregates length
% 16: m         number of iterations
VarNameArray = char('klpf','klmf','klpm','klmm','kap','kam','kn','cnmRatio','kbf','kbm','y00','p00','rs','rfm');

% select number of base parameter
basePnum = 11;
% select adjustable parameter
adjPnum = 14;
% Data file protperties
FileName = 'PhaseLine vary y0 and rfm-rs.txt';
FilePath = 'C:\Users\youval\Documents\Education\Research\Prion line of defence\Paper 1\Matlab workspace\Critical aggregation phase line\run_2014';
% FilePath = 'C:\Users\youval\Documents\Education\Research\Prion line of defence\Paper 1\Matlab workspace\40% a-beta reduction phase line\';
% Use this when changing both rfm and rs  DON'T FORGET TO CHECK InitParameterUpdate to see if rs and rfm are changed simultanuously
% FilePath = 'C:\Users\youval\Documents\Education\Research\Prion line of defence\Paper 1\Matlab workspace\40% a-beta reduction phase line\Test runs\';
% select direction to explor (-1) decrease basePval (1) increase
% ExplorDir is being set at FItFunction3
% Plot parameters
xnamestr = VarNameArray(adjPnum,:);
ynamestr = VarNameArray(basePnum,:);


% Initial Values           
ndataPoints = 25;               % the number of new data points we want to fine on the phase path
JobStatus = 1;                  % Job Errors flag
TestVal = zeros(1,4);           % pass test 0, fail test 1
AccomulationFlag = 0;           % check if we have accomulation at the end of the y or g arrays
y = zeros(1000,1);              % Initializing point arrays
x = zeros(1000,1);              % Initializing point arrays
tempX = 0;
tempY = 0;
OldTransPosSign = 0;            % 
OldXY = 0;
dy = 0;
xFactor = 0;
yFactor = 0;
fittedFun = '';
regularSol = 1;                 % when 1 we fit function x=f(y) when 0 we fit y=f(x)

% Open a file for data collection
% the file need to have at leat one line of data for the initial conditions
% Every line in the file that starts with # will be ignorred
fileID = fopen([FilePath FileName],'r+');       % Open file for reading
disp(['Reading and writing data to ' FileName])
disp(['Working on  ' xnamestr ' and ' ynamestr])
oldDataLength = 0;
DataLine = fgetl(fileID);
while ischar(DataLine)
    if DataLine(1)=='['          % Collecting data only from lines that starts with [
%         oldDataLength = oldDataLength + 1;
%         initParameter = eval(DataLine);
%         y(oldDataLength) = initParameter(basePnum);      % Set initial value
%         x(oldDataLength) = initParameter(adjPnum);       % Set initial value    
        % starting fresh  - ignoring old runs
        oldDataLength = oldDataLength+1;
        initParameter = eval(DataLine);
        y(1) = initParameter(basePnum);      % Set initial value
        x(1) = initParameter(adjPnum);       % Set initial value 
        
    end
%     disp(DataLine)
    DataLine = fgetl(fileID);
end

% Start labs
if matlabpool('size')==0
    matlabpool open 4               % allocate local MATLAB labs for parallel functions execution  
else
    matlabpool('close')
    matlabpool open 4 
end
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

if oldDataLength==0
    % no initial point or points
    JobStatus = 0;
    disp('No initial point or points where found in the data file !!!')
    matlabpool close
    fclose(fileID);
end

if JobStatus
    % Start main loop
    % At this point we need to have the old data, if exist in the y and x arrays
    % And the correct file need to be open for writing
%     i = oldDataLength;
    i = 1;
    Extra_n = 1;                % in use when larger n is needed
    Reslimit = 30*pi/180/5;     % angular resolution of the accuracy of the process  (12, 30)
    UpperResLimit = pi/3;       % Limit how large the search angle can be
    initResFactor = pi/3;
    stdResFactor = 30*pi/180;
    [xFactor,yFactor] = XYfactor(x(1),y(1));
    % plot parameters
    xpowerstr = num2str(floor(log10(xFactor)));
    ypowerstr = num2str(floor(log10(yFactor)));
    strx = [xnamestr ' [\times 10^{ ' xpowerstr '}]'];
    stry = [ynamestr ' [\times 10^{ ' ypowerstr '}]'];
    disp('-----   Main loop start  -----')
    while (i < (1+ndataPoints)) && JobStatus;
        disp('*********************************************************************************')
        tic;
        RunCycle = 1;           % If there are no additional run issues, Two run cycles will give a point on the phase line
        UseFitFlag = 1;         % if use fit flag=1 -> use fit fanction, otherwise skip fitting
        disp(['Step number (i): ' num2str(i)])
        % Start looking for new point
        while RunCycle < 3
            if RunCycle==1
                disp(['RunCycle = 1 , UseFitFlag = ' num2str(UseFitFlag)])
                % Clearing history
                TransPos = 0;           % indicate where is the transition between values 0 and 1
                x1 = 0;y1 = 0;x2 = 0;y2 = 0;
                if UseFitFlag
                    % best curve fit using the available data
                    % Calculate r.  We want to have as large r as possible
                    [nFit,fitParameter,fittedFun,regularSol,ExplorDir,r,baseAngle]=FitFunction4(i,x,y,xFactor,yFactor);
                    UseFitFlag = 0;
                    if nFit>1
                        % initial search, at a large angle 
                        resFactor = stdResFactor;
                        disp('Reset resFactor value to stdResFactor')
                    else
                        % Use even larger angle when only one data point is available
                        resFactor = initResFactor;
                    end
                end
                disp(['ExplorDir = ' num2str(ExplorDir)])
                disp(['r = ' num2str(r)])
                disp(['nFit = ' num2str(nFit)])   
            elseif RunCycle==2
                % reduce resolution 
                resFactor = resFactor/5;        
                % continue to converge on the solution till we get to the
                % desiered resolution
                if resFactor <= Reslimit        
                    RunCycle = 3;
                end  
            end
            disp(['resFactor value = ' num2str(resFactor)])
            % Calculate inputMatrix using the linear fit that was found
            % update and prepare the inputMatrix matrix for the 4 labs
            negativeFlag = 1;
            % Ensuring that there are no negative values
            while negativeFlag 
                inputMatrix(1,:) = UpdateVarArray5(initParameter,y(i),x(i),r,baseAngle,xFactor,yFactor,basePnum,adjPnum,ExplorDir,Extra_n,resFactor,1);
                inputMatrix(2,:) = UpdateVarArray5(initParameter,y(i),x(i),r,baseAngle,xFactor,yFactor,basePnum,adjPnum,ExplorDir,Extra_n,resFactor,2);
                inputMatrix(3,:) = UpdateVarArray5(initParameter,y(i),x(i),r,baseAngle,xFactor,yFactor,basePnum,adjPnum,ExplorDir,Extra_n,resFactor,3);
                inputMatrix(4,:) = UpdateVarArray5(initParameter,y(i),x(i),r,baseAngle,xFactor,yFactor,basePnum,adjPnum,ExplorDir,Extra_n,resFactor,4);  
                negativeFlag = testForNegVal(inputMatrix,basePnum,adjPnum);
                if negativeFlag 
                    if r*yFactor>0.001*(y(i))
                        r = r/2;
                        disp('Negative Values - reducing r size')
                    else
                        disp('!!!!!!!!!!!!!!!!!!!')
                        disp('Negative parameters values.   r is very small.')
                        disp('!!!!!!!!!!!!!!!!!!!')
                        break
                    end
                end
            end 
            % Plot existing data and search points
            h = figure(1);
            hold on
            plot(inputMatrix(:,adjPnum)/xFactor,inputMatrix(:,basePnum)/yFactor,'+');
            plot(x(1:i)/xFactor,y(1:i)/yFactor,'-o');
%             plot([x(1:i);xFactor*feval(fittedFun,(y(i)+dy)/yFactor)],[y(1:i);y(i)+dy],'-o');
            drawnow
            % Saving the plot to a file
            filename = [FilePath FileName(1:end-4) ' ' num2str(i+oldDataLength) '.jpg'];
            xlabel(strx);
            ylabel(stry);
            title(FileName(1:end-4))
%             saveas(h,filename);
            hold off
            close(h)
 
            % test inputMatrix size
            % This program designed to run on 4 labs
            if size(inputMatrix,1)~=4
                disp('ERROR - inputMatrix has wrong dimentions')
                JobStatus = 0;
                break
            end
            % Run the job on 4 labs 
            disp('Start running Labs')
            spmd
                [yLast,gLast,tLast,TestFlag,TestFlagValue,n10,n90,x1test,x2test,alarmFlag,t,ft,F,M] = AggKineticsV16PhaseExplor(inputMatrix(labindex,:));
            end
            
            % Plot aggregation data
            for fn=1:4
                h = figure(fn+1);
                xdata = t{fn};
                ydata = ft{fn};
                xx2 = [xdata(n10{fn}) xdata(n10{fn})];
                xx3 = [xdata(n90{fn}) xdata(n90{fn})];
                yy2 = [0 max(ydata)];
                plot(xdata,ydata,xx2,yy2,xx3,yy2)
                saveas(h,['A-Beta-' num2str(fn) '.jpg']);
%                 plot(xdata,ydata,xx2,yy2,2*xx2,yy2,xx3,yy2)
%                 plot(xdata(n10{fn}:end),ydata(n10{fn}:end))                
%                 plot(F{fn})
                drawnow
            end
            % display test results
            disp('TestFlag')
            disp([num2str(TestFlag{1}) '  :  ' num2str(TestFlagValue{1}) ' :: ' num2str([n10{1},n90{1},x1test{1},x2test{1}])])
            disp([num2str(TestFlag{2}) '  :  ' num2str(TestFlagValue{2}) ' :: ' num2str([n10{2},n90{2},x1test{2},x2test{2}])])
            disp([num2str(TestFlag{3}) '  :  ' num2str(TestFlagValue{3}) ' :: ' num2str([n10{3},n90{3},x1test{3},x2test{3}])])
            disp([num2str(TestFlag{4}) '  :  ' num2str(TestFlagValue{4}) ' :: ' num2str([n10{4},n90{4},x1test{4},x2test{4}])])
            
            % Colect the test results and check that we used a large enough
            % array for our simulation. We don't want aggregates at the end
            % of the arrays x,y
            [TestVal,AccomulationFlag] = GetTestValues(TestVal,TestFlag,yLast,gLast,alarmFlag);
            disp(['TestVal = ' num2str(TestVal)])
            disp(['AccomulationFlag = ' num2str(AccomulationFlag)])
            % Get test transion position and store the old position
            OldTransPosSign = sign(TransPos);
            TransPos = GetTransionPos2(TestVal);
            disp(['TransPos = ' num2str(TransPos)])  
            % Set current solution values
            [x1,y1,x2,y2,angleFactor] = SetTempXY(inputMatrix,TransPos,adjPnum,basePnum,x1,y1,x2,y2,OldTransPosSign);
            % Change the base angle to the center of the solution region
            baseAngle = baseAngle + angleFactor*resFactor;
            
            % test initial run and progres to RunCycle 2
            if RunCycle == 1;
                if abs(TransPos)~=4                 % OK to continue
                    if IsIn(abs(TransPos),[1,3])    % if we are not at the middle -> reduce r
                        r = r/2;
                    end
                    RunCycle = 2;
                else                            % 1st run problems
                    disp('Problem with initial run - Reduce r zise and increase angular search resolution')
                    if resFactor<UpperResLimit
                        disp('Adjusting the resFactor angle')
                        resFactor = 2*resFactor;
                        resFactor = resFactor - (resFactor>UpperResLimit)*(resFactor-UpperResLimit);
                        disp(['resFactor = ' num2str(resFactor*pi/180)])
                    else
                        r = r/2;
                        if r*yFactor<0.001*(y(i))
                            disp('r is smaller than 0.1% of y(i)   ->  Stop program')
                            JobStatus = 0;
                            break
                        end
                    end
                end
            end
            % Check if we have aggregates at the end of arrays x,y
            if AccomulationFlag         
                Extra_n = Extra_n*2;    % increase arrays size
                RunCycle = 1;           % Start from the begining
                disp('xxxxxxxxxxxxxxxxxxxxx')
                disp('Need to run again with loger aggregate lenght')
                disp(['i = ' num2str(i)])
                disp(['Extra_n = ' num2str(Extra_n)])
                disp('xxxxxxxxxxxxxxxxxxxxx')
            end
            
            % store current result and go to next point 
            if RunCycle==3          
                % take the center of the solution region
                x(i+1) = (x1+x2)/2;
                y(i+1) = (y1+y2)/2;
                % set and write new data to file
                initParameter = InitParameterUpdate(initParameter,basePnum,y(i+1),adjPnum,x(i+1));
                fprintf(fileID,'%s\n',mat2str(initParameter));
                disp('============================')
                disp(['Recording point number ' num2str(i+oldDataLength)])
                disp(['y(i+1) = ' num2str(y(i+1))])
                disp(['x(i+1) = ' num2str(x(i+1))])
                % Increse i
                i = i+1;
                % Making sure we have no negative values
               if (y(i)<0 || x(i)<0)
                   JobStatus = 0;
                   disp('We have a negative concentration value -> stop program')
                   fclose(fileID);
                   matlabpool close   
               end
               PoitCalcTime = toc;
               disp(['Cycle time: ' num2str(PoitCalcTime)])
            end
            disp(['RunCycle = ' num2str(RunCycle)])           
        end             % end the while RunCycle<3 statement
    end

    % Close Parallel labs
    if matlabpool('size')~=0
        matlabpool close
    end
    % Close file
    fclose(fileID);   
end

PlotPhasePath2(FileName,FilePath,basePnum,adjPnum)


