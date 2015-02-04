function PlotPhaseCollection()
    % This program is an assembly of several pieces of the phase 
    % exploration files. I put them together, to make it easier for other 
    % users to modify and adjust the way it runs.

    % author Youval Dar
    % collection assembly at 4/8/2014

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
    disp('--- Start ---')
    % change random numbers stream
    RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));

    % Initial time scale settings
    % tDay = 60*60*24;    % Seconds in a day
    % tHour = 3600;       % Seconds in one hour
    % tFactor = 1e-6;     % Time factor for simulation
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
    % 17: y0end     A-Beta monomers concentration at the end
    % 18: F_end     A-Beta concentration in aggregate, at the end
    VarNameArray = char('klpf','klmf','klpm','klmm','kap','kam','kn','cnmRatio','kbf','kbm','y00','p00','rs','rfm','AB_average','F_average');

    % select number of base parameter
    basePnum = 11;
    % select adjustable parameter
    adjPnum = 14;
    % Data file protperties
    FileName = 'PhaseLine vary y0 and rfm.txt';
    FilePath = 'C:\Users\Youval\Documents\Education\Research\Prion line of defence\Matlab programs\For Ben\';
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
    y0end = zeros(1000,1);          % Initializing point arrays
    F_end = zeros(1000,1);
    tempX = 0;
    tempY = 0;
    OldTransPosSign = 0;            % 
    OldXY = 0;
    dy = 0;
    xFactor = 0;
    yFactor = 0;
    fittedFun = '';
    regularSol = 1;    % when 1 we fit function x=f(y) when 0 we fit y=f(x)

    % Open a file for data collection
    % the file need to have at leat one line of data for the initial conditions
    % Every line in the file that starts with # will be ignorred
    fileID = fopen([FilePath FileName],'r+');       % Open file for reading
    disp(['Reading and writing data to ' FileName])
    disp(['Working on  ' xnamestr ' and ' ynamestr])
    oldDataLength = 0;
    DataLine = fgetl(fileID);
    while ischar(DataLine)
        % Collecting data only from lines that starts with [
        if DataLine(1)=='['          
            oldDataLength = oldDataLength + 1;
            initParameter = eval(DataLine);

            % Choose if you want to start a new run or continue a previous one
            % starting fresh  - ignoring old runs
            startPos = 1;
            % starting at the end of the old run
    %         startPos = oldDataLength;

            y(startPos) = initParameter(basePnum);      % Set initial value
            x(startPos) = initParameter(adjPnum);       % Set initial value 
            y0end(startPos) = initParameter(17);
            F_end(startPos) = initParameter(18);

        end
    %     disp(DataLine)
        DataLine = fgetl(fileID);
    end

    % Start labs
    disp ('--- Start labs ---')
    if matlabpool('size')==0
        % allocate local MATLAB labs for parallel functions execution  
        matlabpool open 4   
    else
        % close and reopen if labs are already open
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
        i = startPos;
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
                        [nFit,fitParameter,fittedFun,regularSol,ExplorDir,r,baseAngle]=FitFunction(i,x,y,xFactor,yFactor);
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
                    inputMatrix(1,:) = UpdateVarArray(initParameter,y(i),x(i),y0end(i),F_end(i),r,baseAngle,xFactor,yFactor,basePnum,adjPnum,ExplorDir,Extra_n,resFactor,1);
                    inputMatrix(2,:) = UpdateVarArray(initParameter,y(i),x(i),y0end(i),F_end(i),r,baseAngle,xFactor,yFactor,basePnum,adjPnum,ExplorDir,Extra_n,resFactor,2);
                    inputMatrix(3,:) = UpdateVarArray(initParameter,y(i),x(i),y0end(i),F_end(i),r,baseAngle,xFactor,yFactor,basePnum,adjPnum,ExplorDir,Extra_n,resFactor,3);
                    inputMatrix(4,:) = UpdateVarArray(initParameter,y(i),x(i),y0end(i),F_end(i),r,baseAngle,xFactor,yFactor,basePnum,adjPnum,ExplorDir,Extra_n,resFactor,4);  
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
                     [yLast,gLast,tLast,TestFlag,TestFlagValue,n10,n90,x1test,x2test,alarmFlag,t,ft,AB_average,F_average]= AggKineticsV16PhaseExplorCollection(inputMatrix(labindex,:));
%                    [yLast,gLast,tLast,TestFlag,TestFlagValue,n10,n90,x1test,x2test,alarmFlag,t,ft,AB_average,F_average]= AggKineticsV16PhaseExplorCollection(inputMatrix(1,:));
                end
     
                % Plot aggregation data
                for fn=1:4
                    h = figure(fn+1);
                    xdata = t{fn};
                    ydata = ft{fn};
%                     y_AB_average = AB_average{fn};
%                     F_F_average = F_average{fn};
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
                [x1,y1,x2,y2,angleFactor,ab_average,f_average] = SetTempXY(inputMatrix,TransPos,adjPnum,basePnum,x1,y1,x2,y2,OldTransPosSign,AB_average,F_average);
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
                    y0end(i+1) = ab_average;
                    F_end(i+1) = f_average;
                    fprintf(fileID,'%s\n',mat2str(initParameter));
                    % set and write new data to file
                    initParameter = InitParameterUpdate(initParameter,basePnum,y(i+1),adjPnum,x(i+1),y0end(i+1),F_end(i+1));
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

    PlotPhasePath(FileName,FilePath,basePnum,adjPnum)
end

function [nFit,fitParameter,fittedFun,regularSol,ExplorDir,r,baseAngle]=FitFunction(i,x,y,xFactor,yFactor)

    % Find the best function to predict the next point
    % do not use more that 6 points
    % set default values
    yy = y/yFactor;
    xx = x/xFactor;
    testvalue = 0.8;
    fitParameter = 0;            % adjrsquare
    nFit = 1;
    baseAngle = pi/2;             % Typically pi/2 but can use 3*pi/4,7*pi/8 and so on
    fittedFun = '';
    regularSol = 1;
    ExplorDir = -1;
    minPoints = 5;      % The minimum number of points we want to consider for the linear fit is minPoints + 1

    if i>minPoints
        j=(i-minPoints);
        % Using only linear fit and only with minPoints + 1 points
        [TempFun,gof]=fit(yy(j:i),xx(j:i),'poly1');
        TempfitParameter = abs(GetFitParameter(gof));

        if TempfitParameter > testvalue
            fitParameter = TempfitParameter;
            fittedFun = TempFun;
            nFit = i-j+1;
            % check if the line is close to horizontal, if yes use horizontal solution
            if sum(yy(i)==yy(j:i))~=(nFit)
                regularSol = 1; 
            else                    % we have a stright horizontal line  
                regularSol = 0;
            end
        end
    end
    % Calculate dy
    if ~isempty(fittedFun)              % Check if we had a good fit 
        dy1 = y(i)-y(i-nFit+1);
        ExplorDir = sign(dy1) + (dy1==0);
        dy = getdy(ExplorDir,y(i),dy1,yFactor,fittedFun,regularSol);
        % Calculate baseAngle
        dyy = dy/yFactor;                                               % we evaluated the function with factorized arrays values
        if regularSol
            newX = feval(fittedFun,(y(i)+dy)/yFactor)*xFactor;          % get the x value
            dx = newX - x(i);
            baseAngle = atan(differentiate(fittedFun,y(i)/yFactor));    % the angle of the perpendicullar line at point y
        else                % Horizontal solution
            dx = 0.1*x(i);
            baseAngle = pi/2;  % the angle of the perpendicullar line at point y
        end
        dxx = dx/xFactor;
        r = (0.5+rand())*sqrt(dyy^2+dxx^2);
    else  
        r = 0.05*y(i)/yFactor;
    end

end

function OutArray = UpdateVarArray(initParameter,y,x,AB_average,F_average,r,baseAngle,xFactor,yFactor,basePnum,adjPnum,ExplorDir,Extra_n,resFactor,labNumber)

    % This function updates the values sent to each lab running the
    % AggKineticsV16PhaseExplorer.m
    % this function is used by ExplorPhaseSpace6.m
    % Note that r is in the factorized units

    fArray = [-1.5,-0.5,0.5,1.5]*resFactor;            % Initial search
    dyy = r*cos(baseAngle+fArray(labNumber));
    dxx = r*sin(baseAngle+fArray(labNumber));
    dy = ExplorDir*dyy*yFactor;
    dx = ExplorDir*dxx*xFactor;
    newY = y+dy;
    newX = x+dx;

    OutArray = InitParameterUpdate(initParameter,basePnum,newY,adjPnum,newX,AB_average,F_average);
    OutArray(15) = round(OutArray(15)*Extra_n);        % Adjust n (aggregate length) if needed

end

function [TestVal,AccomulationFlag] = GetTestValues(TestVal,TestFlag,yLast,gLast,alarmFlag)
    % Get test values for PlotPhaseCollection
    % TestVal is an array (1,4). it idicates wether a test passed or not. pass test 0, fail test 1
    % AccomulationFlag check if we have accomulation at the end of the y or x arrays
    % If we have aggregates at the and we need to increase the maximum size of y and x arrays 

    AccomulationFlag = 0;

    for j=1:4
        % test fail if we have aggregates for at leat FailCounter steps
        TestVal(j) = TestFlag{j};             
        AccomulationFlag = AccomulationFlag + (yLast{j}>0 || gLast{j}>0);
        disp(['Lab number: ' num2str(j)])
        disp('==========================')
        disp(['Alarm Value (Particle number factor): ' num2str(alarmFlag{j})])
        disp(['Test Value: ' num2str(TestFlag{j})])
        disp(['yLast: ' num2str(yLast{j})])
        disp(['gLast: ' num2str(gLast{j})])
        disp('--------------------------')
    end
    disp(['AccomulationFlag: ' num2str(max(AccomulationFlag))])
    disp('==========================')
end

function [xFactor,yFactor] = XYfactor(x,y)

    xFactor = 1;
    yFactor = 1;
    % factorizing the x and y so that the fitting function will work better
    if x~=0
        xFactor = 10^floor(log10(x));      
    end
    if y~=0
        yFactor = 10^floor(log10(y));
    end
end

function OutArray = InitParameterUpdate(initParameter,basePnum,NewY,adjPnum,NewX,AB_average,F_average)

    % This function adjusts initParameter array
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
    % 17: y0end     final A-Beta concentration
    % 18: F_end     A-Beta concentration in aggregate, at the end
    % ??: FailCounter

    % Get the previous values
    OutArray = initParameter;
    
    % Update the AB concentration average at the end to the test
    OutArray(17) = AB_average;
    OutArray(18) = F_average;

    % update values
    flag1 = (basePnum == 11 || basePnum == 12);

    % When keeping klmf and klmm equal
    flag2 = (adjPnum == 2 || adjPnum == 4);
    % When keeping kbf and  kbm equal
    flag3 = (adjPnum == 9 || adjPnum == 10);
    % When keeping klpf and klpm equal
    flag4 = (adjPnum == 1 || adjPnum == 3);
    % When keeping rs and rfm equal
    flag5 = (adjPnum == 13 || adjPnum == 14);

    if  flag1 && flag2 
        OutArray(basePnum) = NewY;
        OutArray(2) = NewX;
        OutArray(4) = NewX;
    elseif flag1 && flag3
        OutArray(basePnum) = NewY;
        OutArray(9) = NewX;
        OutArray(10) = NewX;
    elseif flag1 && flag4
        OutArray(basePnum) = NewY;
        OutArray(1) = NewX;
        OutArray(3) = NewX;
    % elseif flag1 && flag5
    %     OutArray(basePnum) = NewY;
    %     OutArray(13) = NewX;
    %     OutArray(14) = NewX;
    else
        OutArray(basePnum) = NewY;
        OutArray(adjPnum) = NewX;
    end
end

function [negativeFlag] = testForNegVal(inputMatrix,basePnum,adjPnum)

    % checking if any variables are negative
    % Parameters numbers in inputMatrix
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

    negativeFlag = 0;

    for i=1:4
        if inputMatrix(1,basePnum)<0 || inputMatrix(1,adjPnum)<0
            negativeFlag = 1;
        end
    end
end

function [x1,y1,x2,y2,angleFactor,ab_average,f_average] = SetTempXY(inputMatrix,TransPos,adjPnum,basePnum,x1,y1,x2,y2,OldTransPosSign,AB_average,F_average)
    % angleFactor indicates how the base angle should be adjusted, direction
    % and units of angle resulotion at the particular run
    switch TransPos         % done with this point - set values of y0,x0
        case -1             % [1,0,0,0] 
            x1 = inputMatrix(1,adjPnum);    % Set current value      
            y1 = inputMatrix(1,basePnum);   % Set current value  
            x2 = inputMatrix(2,adjPnum);    % Set current value      
            y2 = inputMatrix(2,basePnum);   % Set current value
            angleFactor = -1;
            ab_average = AB_average{1};
            f_average = F_average{1};
        case 1              % [0,1,1,1]
            x1 = inputMatrix(1,adjPnum);    % Set current value
            y1 = inputMatrix(1,basePnum);   % Set current value
            x2 = inputMatrix(1,adjPnum);    % Set current value
            y2 = inputMatrix(1,basePnum);   % Set current value
            angleFactor = -1;
            ab_average = AB_average{2};
            f_average = F_average{2};
        case -2             % [1,1,0,0]
            x1 = inputMatrix(2,adjPnum);    % Set current value      
            y1 = inputMatrix(2,basePnum);   % Set current value
            x2 = inputMatrix(3,adjPnum);    % Set current value      
            y2 = inputMatrix(3,basePnum);   % Set current value
            angleFactor = 0;
            ab_average = AB_average{2};
            f_average = F_average{2};
        case 2              % [0,0,1,1]
            x1 = inputMatrix(2,adjPnum);    % Set current value
            y1 = inputMatrix(2,basePnum);   % Set current value
            x2 = inputMatrix(2,adjPnum);    % Set current value
            y2 = inputMatrix(2,basePnum);   % Set current value
            angleFactor = 0;
            ab_average = AB_average{3};
            f_average = F_average{3};
        case -3             % [1,1,1,0]
            x1 = inputMatrix(3,adjPnum);    % Set current value
            y1 = inputMatrix(3,basePnum);   % Set current value
            x2 = inputMatrix(4,adjPnum);    % Set current value
            y2 = inputMatrix(4,basePnum);   % Set current value
            angleFactor = 1;
            ab_average = AB_average{3};
            f_average = F_average{3};
        case 3              % [0,0,0,1]
            x1 = inputMatrix(3,adjPnum);    % Set current value
            y1 = inputMatrix(3,basePnum);   % Set current value
            x2 = inputMatrix(4,adjPnum);    % Set current value
            y2 = inputMatrix(4,basePnum);   % Set current value
            ab_average = AB_average{4};
            f_average = F_average{4};
            angleFactor = 1;
        case -4             % [1,1,1,1]
            if OldTransPosSign >0
                % x1 stay the same
                % y1 stay the same
                x2 = inputMatrix(1,adjPnum);    % Set current value
                y2 = inputMatrix(1,basePnum);   % Set current value
                angleFactor = -2;
                ab_average = AB_average{1};
                f_average = F_average{1};
            elseif OldTransPosSign <0
                x1 = inputMatrix(4,adjPnum);    % Set current value
                y1 = inputMatrix(4,basePnum);   % Set current value
                % x2 stay the same
                % y2 stay the same
                angleFactor = 2;
                ab_average = AB_average{4};
                f_average = F_average{4};
            else
                disp('OldTransPosSign = 0   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                x1 = 0;
                y1 = 0;
                x2 = 0;
                y2 = 0;
                angleFactor = 0;
                ab_average = 0;
                f_average = 0;
            end   
        case 4              % [0,0,0,0]
            if OldTransPosSign >0
                x1 = inputMatrix(4,adjPnum);    % Set current value
                y1 = inputMatrix(4,basePnum);   % Set current value
                % x2 stay the same
                % y2 stay the same 
                angleFactor = 2;
                ab_average = AB_average{4};
                f_average = F_average{4};
            elseif OldTransPosSign <0
                % x1 stay the same
                % y1 stay the same
                x2 = inputMatrix(1,adjPnum);    % Set current value
                y2 = inputMatrix(1,basePnum);   % Set current value
                angleFactor = -2;
                ab_average = AB_average{4};
                f_average = F_average{4};
            else
                disp('OldTransPosSign = 0   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                x1 = 0;
                y1 = 0;
                x2 = 0;
                y2 = 0;
                angleFactor = 0;
                ab_average = 0;
                f_average = 0;
            end
    end 
end

function PlotPhasePath(FileName,FilePath,basePnum,adjPnum)
    % % Plot phase path
    % FileName = 'PhaseLine vary y0 and kam.txt';
    % % FilePath = 'C:\Users\youval\Documents\Education\Research\Prion line of defence\Paper 1\Matlab workspace\Critical aggregation phase line\';
    % FilePath = 'C:\Users\youval\Documents\Education\Research\Prion line of defence\Paper 1\Matlab workspace\40% a-beta reduction phase line\';
    % % select number of base parameter
    % basePnum = 11;
    % % select adjustable parameter
    % adjPnum = 6;

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
    % 17: y0end     A-Beta monomer concentration at the end
    % 18: F_end     A-Beta concentration in aggregate, at the end

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
    ParametersName = char('klpf','klmf','klpm','klmm','kap','kam','kn','cnmRatio','kbf','kbm','y00','p00','rs','rfm');



    y = zeros(1000,1);              % Initializing point arrays
    x = zeros(1000,1);              % Initializing point arrays
    y0end = zeros(1000,1);          % Initializing point arrays
    F_end = zeros(1000,1);          % Initializing point arrays

    DataLine = fgetl(fileID);
    i = 1;
    while ischar(DataLine)
        if DataLine(1)=='['          % Collecting data only from lines that starts with [
            ParameterArray = eval(DataLine);
            y(i) = ParameterArray(basePnum);      % Set initial value
            x(i) = ParameterArray(adjPnum);       % Set initial value   
            y0end(i) = ParameterArray(17);
            F_end(i) = ParameterArray(18);
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
    y0end = y0end(1:Endi);
    F_end = F_end(1:Endi);
    % sorting the points
    xy  = sortrows([x y y0end F_end]);
    x = xy(:,1);
    y = xy(:,2);
    y0end = xy(:,3);
    F_end = xy(:,4);

    [xFactor,yFactor] = XYfactor(x(1),y(1));
    x = x(1:Endi)/xFactor;
    y = y(1:Endi)/yFactor;
    y0end = y0end(1:Endi)/yFactor;
    F_end = F_end(1:Endi)/yFactor;
    [fittedFun,fitParameter] = FitFunctionDataPlot(x,y);
    h = figure(1);
    plot(x,y,'o',x,y0end,'o',x,F_end,'o')
    hold on
    plot(fittedFun)
    %plot label
    % xlabel([ParametersName(adjPnum) ' x10^{ ' num2str(floor(log10(xFactor))) '}']);
    % ylabel([ParametersName(basePnum) ' x10^{ ' num2str(floor(log10(yFactor))) '}']);
    xpowerstr = num2str(floor(log10(xFactor)));
    ypowerstr = num2str(floor(log10(yFactor)));
    xnamestr = ParametersName(adjPnum,:);
    ynamestr = ParametersName(basePnum,:);
    strx = [xnamestr ' [\times 10^{ ' xpowerstr '}]'];
    stry = [ynamestr ' [\times 10^{ ' ypowerstr '}]'];
    xlabel(strx);
    ylabel(stry);
    title(FileName(1:end-4))
    if abs(y-1)<0.1
        ylim([0.5 1.5])
    end
    % Saving the plot to a file
    filename = [FilePath FileName(1:end-4) ' data fit  ' '.jpg'];
    saveas(h,filename);
    hold off

    disp fittedFun
end

function [yLast,gLast,tLast,TestFlag,TestFlagValue,n10,n90,x1,x2,alarm,tt,ft,AB_average,F_average]=AggKineticsV16PhaseExplorCollection(DataArray)
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
    % General initialization code for the aggregation functions
    % Use by eval('AggKineticsInitializationV16')

    % Convert to Gillespie and Normalize parameters
    %(klpf,klmf,klpm,klmm,kap,kam,kn,knm,kbf,kbm,y0,p0,rs,n,m,v0)
    v0 = 1e15;          % in nm^3
    na = 6.022e23;      % Avogadro Num
    v1 = v0*1e-24;      % Convert volume from (nm)^3 to litters 
    v = v1*na;          % facilitate conversion of # of mols to # of particles

    clpf = klpf/v;
    clmf = klmf;
    clpm = klpm/v;
    clmm = klmm;
    cap = kap/v;
    cam = kam;
    cn = 2*kn/v;        % free nucleation
    cbf = kbf;          % Free aggregates breaking rate
    cbm = kbm;          % Membrane aggregates breaking rate
    y0 = round(v*y00);  % Convert concentration to number of particles in v0
    p0 = round(v*p00);  % Convert concentration to number of particles in v0
    % rp = rs*p00;
    % non-Gillespie parameters
    crs = rs;       % Prion Sequestration rate
    crp = crs*p0;   % Production rate        
    crfm = rfm;     % Free aggregates out flow rate
    crfp = crfm*y0; % flow supply of free monomers (crfp = crfm*y0)v^2*rfm*y00 = v*rfm*y0


    % Scalling  factor. 
    vScale = 1;
    clpf = clpf*vScale;         % a2 = clpf*y(1)*y(2:(n-1));
    % clmf Stay the same        % a4 = 2*clmf*y(2:n);
    clpm = clpm*vScale;         % a6 = clpm*y(1)*g(1:(n-1));
    % clmm = Stay the same      % a8 = clmm*g(2:n);
    cap = cap*vScale;           % a5 = cap*g0*y;
    % cam Stay the same         % a9 = cam*g(2:n)
    % cn : update at the calculation        % a1 = cn*y(1)*(vScale*y(1)-1)/2;  
    % cbf Stay the same                     % a3 = cbf*y(4:n).*(i-3);     
    % cbm Stay the same                     % a7 = cbm*g(4:n).*(i-2)    
    y0 = round(y0/vScale);      
    p0 = round(p0/vScale);     
    % crs Stay the same
    crp = crp/vScale;           
    % crfm = Stay the same
    crfp = crfp/vScale;
    % t = t*v/3600;   Convert to real time

    if p0~=0
        cnm = vScale*cnmRatio*cn/p0;        % Note the scaling
        % a10 = vScale*cnm*g0*y(1)*(vScale*y(1)-1)/2
    else
        cnm = 0;
    end


    % Initialize parameters

    t = zeros(1,m);         % collect the time steps
    F = zeros(1,m);         % Collect the amount of monomers in the free aggregates at each time step
    M = zeros(1,m);         % Collect the amount of monomers in the membrane aggregates at each time step
    y0t = zeros(1,m);       % the number of free monomers
    g0t = zeros(1,m);       % the number of membrane single proteins
    inagg = zeros(1,m);     % Collect data on sequestered aggregates
    y = zeros(1,n);         % The number of free aggregates of different length
    g = zeros(1,n);         % The number of membrane aggregates of different lengt
    dyin = 0;               % initialize parameters for the deterministic process
    dg0in = 0;              % initialize parameters for the deterministic process
    aa = zeros(1,12);       % Collect the propensity of each region
    % dg0 Collect the fractions of the g0
    % dy1 Collect the fractions of the y1
    % the aggregate clearance is being modified to represent an effective
    % aggregates clearance rate
    crfmArray = crfm*ones(1,n);
    % crfmArray(2:n) = crfmArray(2:n)/50;  
    crsArray = crs*ones(1,n);
    % crsArray(2:n) = crsArray(2:n)/50; 

    % Critical nucleation length in this simulation is 2
    y(1) = y0;              % initial amount of monomers
    g0 = p0;                % initial amount of monomers 
    y0t(1)=y(1);            % tracking free monomers as a function of time
    g0t(1)=g0;              % tracking membrane bound monomers as a function of time
    rst = [1 2 4 2 2 2 4 2 2 1 2 2];            % the start i in each reaction region
    re = [1 (n-1) n n n (n-1) n n n 1 n n];     % the end i in each reaction region


    % initial values for the different reaction propensities
    % Do not include in/out-flow reactions in the Gillespie process  
    % a(1),a(10) are length 1
    a = zeros(12,n);

    % None zero initial values
    a(1,1) = cn*y(1)*(vScale*y(1)-1)/2;
    a(10,1) = cnm*g0*y(1)*(vScale*y(1)-1)/2;

    % Sum all the none zero values
    aa(1) = a(1,1);
    aa(10) = a(10,1);
    a0 = sum(aa);
    % main loop
    % m number of iterations
    i = 2;
    isteps = 1;
    while i<=m && isteps<1e9
        % Add time step
        if a0 ~= 0
            dt = log(1/rand)/a0;
        else
            dt = 1;
        end

        % Get the reaction type number                           
        r = rand;           % Random number (0,1)
        rn = 1;             % the reaction number in a (1 to 12)
        s = aa(rn)/a0;      % accomulating the propensities
        while s < r
            rn = rn + 1;
            s = s + aa(rn)/a0;  
        end
        % find the reaction within the reaction type
        % agn1 is the possition in the y or g arrays
        agn1 = rst(rn);        % the first position in the correct region of the propencity array
        if rn~=1 && rn~=10     % for rn=1 and rn=10 we only have one value in the propencity array                        
            s = s - aa(rn)/a0 + a(rn,agn1)/a0;      % return value to begining of reaction
            while s < r   
                agn1 =agn1 + 1;
                s = s + a(rn,agn1)/a0; 
            end
        end

        % The deterministic flow
        % Adjust y(1) and g0 
        if crfm>0
            dy1 = (crfp/crfm*(1-exp(-crfm*dt)) + y(1)*exp(-crfm*dt)) + dyin - y(1);   % for y(1)
        else
            dy1 = 0;
        end
        dy1Int = fix(dy1);                                                      % Round toward zero
        dyin = dy1 - dy1Int;                                                    % Collecting fractions of particles
        y(1) = y(1) + dy1Int;                                                  % Changing the number of particles 
        y(1) = y(1)*(y(1)>0);                                                % No negative monomer numbers

        % Membrane protein production
        if crs>0
            dg0 = (crp/crs*(1-exp(-crs*dt)) + g0*exp(-crs*dt)) + dg0in - g0;  % for g0
        else
            dg0 = 0;
        end
        dg0Int = fix(dg0);                                              % Round toward zero
        dg0in = dg0 - dg0Int;                                           % Collecting fractions of particles
        g0 = g0 + dg0Int;                                               % Changing the number of particles 
        g0 = g0*(g0>0);                                                 % No negative monomer numbers


        % Calculate the changes according to the reaction that took place   
        % F is the amount of ABeta in free aggregates
        % M is the amount of ABeta in membrane aggregates
        switch rn
            case 1
                % Calculate y and g
                if y(1)>2
                    y(1) = y(1)-2;
                    y(2) = y(2)+1;
                    F(i) = F(i-1)+ 2;
                    % reaction propensities that need update
                    % a1 = cn*y(1)*(y(1)-1)/2;  Updated at the end
                    % a2 = clpf*y(1)*y(2:(n-1));    Updated at the end 
                    % a4 = 2*clmf*y(2:n);
                    a(4,2) = clmf*y(2);
                    aa(4) = sum(a(4,:));
                    % a5 = cap*g0*y;    Updated at the end
                    % a6 = clpm*y(1)*g(1:(n-1));      Updated at the end
                    % a10 = cnm*g0*y(1)*(y(1)-1)/2      Updated at the end
                    % a12 = crfmArray(2:n)*y(2:n);
                    a(12,2) = crfmArray(2)*y(2);
                    aa(12) = sum(a(12,:));
                else
                    F(i) = F(i-1);
                end
                M(i) = M(i-1); 
            case 2
                if y(1)>0 && y(agn1)>0
                    % agn2 is the possition in the y or g arrays
                    agn2 = agn1+1;
                    % Calculate y and g
                    y(agn1) = y(agn1)-1;
                    y(1) = y(1)-1;
                    y(agn2) = y(agn2)+1;
                    F(i) = F(i-1)+ 1;
                    % Calculate changes
                    if agn1>3
                        % a3 = cbf*y(4:n).*(i-3);
                        a(3,agn1) = cbf*y(agn1)*(agn1-3);
                        a(3,agn2) = cbf*y(agn2)*(agn2-3);
                        aa(3) = sum(a(3,:));
                    elseif agn2>3
                        % a3 = cbf*y(4:n).*(i-3);
                        a(3,agn2) = cbf*y(agn2)*(agn2-3);
                        aa(3) = sum(a(3,:));
                    end
                    % reaction propensities that need update
                    % a1 = cn*y(1)*(y(1)-1)/2;     Updated at the end
                    % a2 = clpf*y(1)*y(2:(n-1));    Updated at the end 
                    % a4 = 2*clmf*y(2:n);
                    a(4,agn1) = (1+(agn1>2))*clmf*y(agn1);
                    a(4,agn2) = 2*clmf*y(agn2);
                    aa(4) = sum(a(4,:));
                    % a5 = cap*g0*y;    Updated at the end
                    % a6 = clpm*y(1)*g(1:(n-1));      Updated at the end
                    % a10 = cnm*g0*y(1)*(y(1)-1)/2      Updated at the end
                    % a12 = crfmArray(2:n)*y(2:n);
                    a(12,agn1) = crfmArray(agn1)*y(agn1);
                    a(12,agn2) = crfmArray(agn1)*y(agn2);
                    aa(12) = sum(a(12,:));
                else
                    F(i) = F(i-1);
                end
                M(i) = M(i-1);  
            case 3
                if y(agn1)>0
                    % agn2,3 are the possition in the y or g arrays
                    agn2 = randi(agn1-3)+1;       % The aggregate can break at any bond, randomly
                    agn3 = agn1-agn2;
                    % Calculate y ang g 
                    y(agn1) = y(agn1)-1;               
                    y(agn2) = y(agn2)+1;        % y(ii)-> y(jj)+y(11-jj)
                    y(agn3) = y(agn3)+1;
                    % Calculate changes
                    if agn2~=agn3
                        % a2 = clpf*y(1)*y(2:(n-1));    Updated at the end 
                        % a4 = 2*clmf*y(2:n);
                        a(4,agn3) = (1+(agn3>2))*clmf*y(agn3);
                        aa(4) = sum(a(4,:));
                        % a5 = cap*g0*y;    Updated at the end
                        if agn3>3
                            % a3 = cbf*y(4:n).*yf;
                            a(3,agn3) = cbf*y(agn3)*(agn3-3);
                            aa(3) = sum(a(3,:));
                        end
                    end 
                    if agn2>3
                        % a3 = cbf*y(4:n).*yf;
                        a(3,agn2) = cbf*y(agn2)*(agn2-3);
                        aa(3) = sum(a(3,:));
                    end          
                    % reaction propensities that need update
                    % a2 = clpf*y(1)*y(2:(n-1));    Updated at the end 
                    % a3 = cbf*y(4:n).*yf; 
                    a(3,agn1) = cbf*y(agn1)*(agn1-3);
                    aa(3) = sum(a(3,:));
                    % a4 = 2*clmf*y(2:n);
                    a(4,agn1) = 2*clmf*y(agn1);
                    a(4,agn2) = (1+(agn2>2))*clmf*y(agn2);
                    aa(4) = sum(a(4,:));
                    % a5 = cap*g0*y;    Updated at the end
                    % a12 = crfmArray(2:n)*y(2:n);
                    a(12,agn1) = crfmArray(agn1)*y(agn1);
                    a(12,agn2) = crfmArray(agn2)*y(agn2);
                    a(12,agn3) = crfmArray(agn3)*y(agn3);
                    aa(12) = sum(a(12,:));
                end
                % F(i) does not change
                F(i) = F(i-1);
                M(i) = M(i-1);       
            case 4
                if y(agn1)>0
                    % agn2 is the possition in the y or g arrays
                    agn2 = agn1 - 1;
                    % Calculate y and g
                    y(agn1) = y(agn1)- 1;   
                    y(1) = y(1)+ 1;
                    y(agn2) = y(agn2)+ 1; 
                    if agn1>2
                        F(i) = F(i-1)- 1;
                    else
                        F(i) = F(i-1)- 2;
                    end
                    % Calculate changes
                    if agn1>3
                        % a3 = cbf*y(4:n).*yf;
                        a(3,agn1) = cbf*y(agn1)*(agn1-3);
                        aa(3) = sum(a(3,:));
                    end
                    if agn2>3
                        % a3 = cbf*y(4:n).*yf;
                        a(3,agn2) = cbf*y(agn2)*(agn2-3);
                        aa(3) = sum(a(3,:));
                    end
                    if agn2>1
                        % a4 = 2*clmf*y(2:n);
                        a(4,agn2) = (1+(agn2>2))*clmf*y(agn2);
                        aa(4) = sum(a(4,:));
                        % a5 = cap*g0*y;    Updated at the end
                    end
                    % reaction propensities that need update 
                    % a4 = 2*clmf*y(2:n);
                    a(4,agn1) = (1+(agn1>2))*clmf*y(agn1);
                    aa(4) = sum(a(4,:));
                    % a1 = cn*y(1)*(y(1)-1)/2;   Updated at the end
                    % a2 = clpf*y(1)*y(2:(n-1));    Updated at the end 
                    % a5 = cap*g0*y;    Updated at the end
                    % a6 = clpm*y(1)*g(1:(n-1));      Updated at the end
                    % a10 = cnm*g0*y(1)*(y(1)-1)/2      Updated at the end
                    % a12 = crfmArray(2:n)*y(2:n);
                    a(12,agn1) = crfmArray(agn1)*y(agn1);
                    if agn2>1
                        a(12,agn2) = crfmArray(agn2)*y(agn2);
                    end
                    aa(12) = sum(a(12,:));
                else
                    F(i) = F(i-1);
                end
                M(i) = M(i-1);    
            case 5
                if y(agn1)>0 && g0>0
                    %Adjust y and g
                    y(agn1) = y(agn1)- 1;
                    g0 = g0 - 1;
                    g(agn1) = g(agn1)+ 1;              
                    F(i) = F(i-1)-agn1;  % agn1 should be always >1
                    M(i) = M(i-1)+agn1; 
                    % Calculate Changes         
                    if agn1>3
                        % a3 = cbf*y(4:n).*yf;
                        a(3,agn1) = cbf*y(agn1)*(agn1-3);
                        aa(3) = sum(a(3,:));
                        % a7 = cbm*g(4:n).*(i-2)
                        a(7,agn1)= cbm*g(agn1)*(agn1-2);
                        aa(7) = sum(a(7,:));
                    end
                    % reaction propensities that need update
                    % a2 = clpf*y(1)*y(2:(n-1));    Updated at the end 
                    % a4 = 2*clmf*y(2:n);
                    a(4,agn1) = (1+(agn1>2))*clmf*y(agn1);
                    aa(4) = sum(a(4,:));
                    % a5 = cap*g0*y;    Updated at the end
                    % a6 = clpm*y(1)*g(1:(n-1));      Updated at the end
                    % a8 = clmm*g(2:n);
                    a(8,agn1)= clmm*g(agn1);
                    aa(8) = sum(a(8,:));
                    % a9 = cam*g(2:n);
                    a(9,agn1)= cam*g(agn1);
                    aa(9) = sum(a(9,:));
                    % a10 = cnm*g0*y(1)*(y(1)-1)/2      Updated at the end
                    % a11 = crs*g(1:n)
                    a(11,agn1)= crsArray(agn1)*g(agn1);
                    aa(11) = sum(a(11,:));
                    % a12 = crfmArray(2:n)*y(2:n);
                    a(12,agn1) = crfmArray(agn1)*y(agn1);
                    aa(12) = sum(a(12,:));
                else
                    F(i) = F(i-1);  
                    M(i) = M(i-1);
                end  
            case 6
                if g(agn1)>0 && y(1)>0
                    % agn2 is the possition in the y or g arrays
                    agn2 = agn1 + 1;
                    % Calculate y and g
                    g(agn1) = g(agn1)- 1;
                    y(1) = y(1)- 1;
                    g(agn2) = g(agn2)+ 1;
                    M(i) = M(i-1)+ 1;   
                    % Calculate changes
                    if agn1>3
                        % a7 = cbm*g(4:n).*(i-2)
                        a(7,agn1)= cbm*g(agn1)*(agn1-2);
                        a(7,agn2)= cbm*g(agn2)*(agn2-2);
                        aa(7) = sum(a(7,:));
                    elseif agn2>3
                        % a7 = cbm*g(4:n).*(i-2)
                        a(7,agn2)= cbm*g(agn2)*(agn2-2);
                        aa(7) = sum(a(7,:));
                    end
                    if agn1>1
                        % a8 = clmm*g(2:n);
                        a(8,agn1)= clmm*g(agn1);
                        aa(8) = sum(a(8,:)); 
                    end
                    % reaction propensities that need update
                    % a1 = cn*y(1)*(y(1)-1)/2;     Updated at the end
                    % a2 = clpf*y(1)*y(2:(n-1));    Updated at the end 
                    % a6 = clpm*y(1)*g(1:(n-1));      Updated at the end
                    % a8 = clmm*g(2:n);  agn2>1
                    a(8,agn2)= clmm*g(agn2);
                    aa(8) = sum(a(8,:));
                    % a9 = cam*g(2:n);
                    a(9,agn1)= cam*g(agn1);
                    a(9,agn2)= cam*g(agn2);
                    aa(9) = sum(a(9,:));
                    % a10 = cnm*g0*y(1)*(y(1)-1)/2      Updated at the end
                    % a11 = crs*g(1:n)
                    a(11,agn1)= crsArray(agn1)*g(agn1);
                    a(11,agn2)= crsArray(agn2)*g(agn2);
                    aa(11) = sum(a(11,:));
                else
                    M(i) = M(i-1);
                end
                F(i) = F(i-1);
            case 7
                if g(agn1)>0
                    % agn2,3 are the possition in the y or g arrays
                    agn2 = randi(agn1-3)+1;       % The aggregate can break at any bond, randomly (g length)
                    agn3 = agn1-agn2;           % ang3 relate to the y length
                    % Calculate y and g
                    g(agn1) = g(agn1)- 1; 
                    g(agn2) = g(agn2)+ 1;
                    y(agn3) = y(agn3)+ 1;  
                    F(i) = F(i-1)+agn3;
                    M(i) = M(i-1)-agn3; 
                    % Calculate changes
                    if agn3>3
                        % a3 = cbf*y(4:n).*yf;
                        a(3,agn3) = cbf*y(agn3)*(agn3-3);
                        aa(3) = sum(a(3,:));
                    end
                    if agn3>1
                        % a4 = 2*clmf*y(2:n);
                        a(4,agn3) = (1+(agn3>2))*clmf*y(agn3);
                        aa(4) = sum(a(4,:));
                    end
                    if agn2>1
                        % a8 = clmm*g(2:n);
                        a(8,agn2)= clmm*g(agn2);
                        aa(8) = sum(a(8,:));
                    end
                    if agn2>3
                        % a7 = cbm*g(4:n).*(i-2)
                        a(7,agn2)= cbm*g(agn2)*(agn2-2);
                        aa(7) = sum(a(7,:));
                    end  
                    % reaction propensities that need update
                    % a2 = clpf*y(1)*y(2:(n-1));    Updated at the end  
                    % a5 = cap*g0*y;    Updated at the end
                    % a6 = clpm*y(1)*g(1:(n-1));      Updated at the end
                    % a7 = cbm*g(4:n).*(i-2)
                    a(7,agn1)= cbm*g(agn1)*(agn1-2);
                    aa(7) = sum(a(7,:));
                    % a8 = clmm*g(2:n);
                    % a8 = clmm*g(2:n);
                    a(8,agn1)= clmm*g(agn1);
                    aa(8) = sum(a(8,:));
                    % a9 = cam*g(2:n);
                    a(9,agn1)= cam*g(agn1);
                    a(9,agn2)= cam*g(agn2);
                    aa(9) = sum(a(9,:));
                    % a11 = crs*g(1:n)
                    a(11,agn1)= crsArray(agn1)*g(agn1);
                    a(11,agn2)= crsArray(agn2)*g(agn2);
                    aa(11) = sum(a(11,:));
                    % a12 = crfmArray(2:n)*y(2:n);
                    if agn3>1
                        a(12,agn3) = crfmArray(agn3)*y(agn3);
                    end
                    aa(12) = sum(a(12,:));
                else
                    F(i) = F(i-1);
                    M(i) = M(i-1);
                end 
            case 8
                if g(agn1)>0
                    % agn2 is the possition in the y or g arrays
                    agn2 = agn1 - 1; 
                    % Calculate y and g
                    if agn1>2
                        g(agn1) = g(agn1)- 1;   
                        g(agn2) = g(agn2)+ 1; 
                        y(1) = y(1)+ 1;
                        M(i) = M(i-1)- 1;   
                    else
                        % if agn1=2 it falls apart to 2 seperate monomers
                        g(agn1) = g(agn1)- 1;   
                        g0 = g0+ 1; 
                        y(1) = y(1)+ 2;
                        M(i) = M(i-1)- 2; 
                    end
                    % Calculate changes
                    if agn1>3
                        % a7 = cbm*g(4:n).*(i-2)
                        a(7,agn1)= cbm*g(agn1)*(agn1-2);
                        aa(7) = sum(a(7,:));
                    end
                    if agn2>3
                        % a7 = cbm*g(4:n).*(i-2)
                        a(7,agn2)= cbm*g(agn2)*(agn2-2);
                        aa(7) = sum(a(7,:));
                    end
                    if agn2>1
                        % a8 = clmm*g(2:n);
                        a(8,agn2)= clmm*g(agn2);
                        aa(8) = sum(a(8,:));
                    end
                    % reaction propensities that need update
                    % a1 = cn*y(1)*(y(1)-1)/2;     Updated at the end
                    % a2 = clpf*y(1)*y(2:(n-1));    Updated at the end    
                    % a6 = clpm*y(1)*g(1:(n-1));      Updated at the end  
                    % a8 = clmm*g(2:n);
                    a(8,agn1)= clmm*g(agn1);
                    aa(8) = sum(a(8,:));
                    % a9 = cam*g(2:n);
                    a(9,agn1)= cam*g(agn1);
                    if agn2>1
                        a(9,agn2)= cam*g(agn2);
                    end
                    aa(9) = sum(a(9,:));
                    % a10 = cnm*g0*y(1)*(y(1)-1)/2      Updated at the end
                    % a11 = crs*g(1:n)
                    a(11,agn1)= crsArray(agn1)*g(agn1);
                    if agn2>1
                        a(11,agn2)= crsArray(agn2)*g(agn2);
                    end
                    aa(11) = sum(a(11,:));
                else
                    M(i) = M(i-1);
                end
                F(i) = F(i-1);
            case 9
                if  g(agn1)>0
                    % Calculate y and g
                    y(agn1) = y(agn1)+ 1;
                    g0 = g0 + 1;
                    g(agn1) = g(agn1)- 1;
                    if agn1>1
                        F(i) = F(i-1)+ agn1;
                    else
                        F(i) = F(i-1);
                    end
                    M(i) = M(i-1)- agn1;
                    % Calculate changes
                    % a2 = clpf*y(1)*y(2:(n-1));    Updated at the end 
                    % a4 = 2*clmf*y(2:n);
                    a(4,agn1) = (1+(agn1>2))*clmf*y(agn1);
                    aa(4) = sum(a(4,:));
                    % a6 = clpm*y(1)*g(1:(n-1));      Updated at the end
                    % a8 = clmm*g(2:n);
                    a(8,agn1)= clmm*g(agn1);
                    aa(8) = sum(a(8,:));
                    if agn1>3
                        % a3 = cbf*y(4:n).*yf;
                        a(3,agn1) = cbf*y(agn1)*(agn1-3);
                        aa(3) = sum(a(3,:));
                        % a7 = cbm*g(4:n).*(i-2)
                        a(7,agn1)= cbm*g(agn1)*(agn1-2);
                        aa(7) = sum(a(7,:));
                    end         
                    % reaction propensities that need update
                    % a1 = cn*y(1)*(y(1)-1)/2;     Updated at the end
                    % a5 = cap*g0*y;    Updated at the end
                    % a9 = cam*g(2:n);
                    a(9,agn1)= cam*g(agn1);
                    aa(9) = sum(a(9,:));
                    % a10 = cnm*g0*y(1)*(y(1)-1)/2      Updated at the end
                    % a11 = crs*g(1:n)
                    a(11,agn1)= crsArray(agn1)*g(agn1);
                    aa(11) = sum(a(11,:));
                    % a12 = crfmArray(2:n)*y(2:n);
                    a(12,agn1) = crfmArray(agn1)*y(agn1);
                    aa(12) = sum(a(12,:));
                else
                    F(i) = F(i-1);
                    M(i) = M(i-1);
                end
            case 10
                if y(1)>0 && g0>0 && cnm>0
                    %Adjust y and g
                    y(1) = y(1) - 2;
                    g0 = g0 - 1;
                    g(2) = g(2)+ 1;     
                    M(i) = M(i-1) + 2;   
                    % Calculate changes
                    % reaction propensities that need update
                    % a1 = cn*y(1)*(y(1)-1)/2;      Updated at the end
                    % a2 = clpf*y(1)*y(2:(n-1));    Updated at the end 
                    % a5 = cap*g0*y;    Updated at the end
                    % a6 = clpm*y(1)*g(1:(n-1));      Updated at the end
                    % a8 = clmm*g(2:n);
                    a(8,2)= clmm*g(2);
                    aa(8) = sum(a(8,:));
                    % a9 = cam*g(2:n);
                    a(9,2)= cam*g(2);
                    aa(9) = sum(a(9,:));
                    % a11 = crs*g(1:n)
                    a(11,2)= crsArray(2)*g(2);
                    aa(11) = sum(a(11,:));
                else
                    M(i) = M(i-1);
                end
                F(i) = F(i-1);      % monomers are not counted in F
            case 11
                if g(agn1)>0
                    g(agn1) = g(agn1) - 1;
                    M(i-1) = M(i-1) - agn1;
                    % Do not do anything to F since i did not change
                    % a6 = clpm*y(1)*g(1:(n-1));      Updated at the end
                    if agn1>3
                        % a7 = cbm*g(4:n).*(i-2)
                        a(7,agn1)= cbm*g(agn1)*(agn1-2);
                        aa(7) = sum(a(7,:));
                    end  
                    % a8 = clmm*g(2:n);
                    a(8,agn1)= clmm*g(agn1);
                    aa(8) = sum(a(8,:));
                    % a9 = cam*g(2:n);
                    a(9,agn1)= cam*g(agn1);
                    aa(9) = sum(a(9,:));
                    % a11 = crs*g(1:n)
                    a(11,agn1)= crsArray(agn1)*g(agn1);
                    aa(11) = sum(a(11,:));
                end
            case 12
                if y(agn1)>0
                    y(agn1) = y(agn1) - 1;
                    F(i-1) = F(i-1) - agn1;
                    % Do not do anything to M since i did not change
                    % Calculate changes
                    % reaction propensities that need update
                    % a1 = cn*y(1)*(y(1)-1)/2;      Updated at the end
                    % a2 = clpf*y(1)*y(2:(n-1));    Updated at the end
                    if agn1>3
                        % a3 = cbf*y(4:n).*(i-3);
                        a(3,agn1) = cbf*y(agn1)*(agn1-3);
                        aa(3) = sum(a(3,:));
                    end
                    % a4 = 2*clmf*y(2:n);
                    a(4,agn1) = (1+(agn1>2))*clmf*y(agn1);
                    aa(4) = sum(a(4,:));
                    % a5 = cap*g0*y;    Updated at the end
                    % a12 = crfmArray(2:n)*y(2:n);
                    a(12,agn1) = crfmArray(agn1)*y(agn1);
                    aa(12) = sum(a(12,:));  
                end    
            otherwise
                % disp('Problem with the reaction number')   
        end

        if a(11,1)~=0
            disp('a(11,1)~=0')
            break
        end

        % Always update 
        y0t(i) = y(1);
        g0t(i) = g0;
        % reaction propensities that need update
        % a1 = cn*y(1)*(y(1)-1)/2;  
        a(1,1) = cn*y(1)*(vScale*y(1)-1)/2;
        aa(1) = a(1,1);
        % a2 = clpf*y(1)*y(2:(n-1));
        a(2,rst(2):re(2))= clpf*y(1)*y(2:(n-1));
        aa(2) = sum(a(2,:));
        % a5 = cap*g0*y;
        % Change the value of cap to favor attachment of particular aggregate length  
        a(5,rst(5):re(5))= cap*g0*y(rst(5):re(5));
    %     AggChangeLength = 10;
    %     a(5,rst(5):AggChangeLength)= cap*g0*y(rst(5):AggChangeLength);
    %     a(5,AggChangeLength+1:re(5))= cap*g0*y(AggChangeLength+1:re(5)).*(AggChangeLength+1:re(5));
    %     a(5,AggChangeLength+1:re(5))= cap*g0*y(AggChangeLength+1:re(5))./(AggChangeLength+1:re(5));
        aa(5) = sum(a(5,:));
        % a6 = clpm*y(1)*g(1:(n-1));
        a(6,rst(6):re(6))= clpm*y(1)*g(2:(n-1));
        aa(6) = sum(a(6,:));
        % a10 = cnm*g0*y(1)*(y(1)-1)/2
        a(10,1) = cnm*g0*y(1)*(vScale*y(1)-1)/2;
        aa(10) = a(10,1); 
        % Sum propensities
        a0 = sum(aa);

        % adjust time steps
        if rn<11
            % Aggregation step
            t(i) = t(i-1)+ dt; 
            i = i + 1;
            isteps = isteps + 1;
        else
            % clearance step
            % do not increase the i counter.
            t(i-1) = t(i-1)+ dt;
            isteps = isteps + 1;
        end
    end

    % convert concentration
    F = F*vScale/v;
    M = M*vScale/v;
    y0t = y0t*vScale/v;
    g0t = g0t*vScale/v;

    % Real world time: treal = t(m)*v/3600;
    % Covert Gillespie time to real time 
    t = t/3600/24;      
    treal = t(end);
    xFactor = floor(log10(t(end)));
    t = t/(10^xFactor);

    if y00<0 || p00<0
        % Test fails -> TestFlag=1
        TestFlag = 1;                   
        TestFlagValue = 1;
    else
        % Test how many of the last 500 steps had segnificant aggregation    
        %  TestVal = 100*F(end-500:end)/(F(end-500:end)+y0t(end-500:end));  

        % for reduced A-Beta use 
    %     [TestFlag,TestFlagValue,n10,n90,x1,x2] = TestABetaValue(y0t,t);

        % For critical aggregation use
        [TestFlag,TestFlagValue,n10,n90,x1,x2,AB_average,F_average] = TestAggValue(F,M,t,y0t);
        %  I take out the last term since in high flow rates it can cause
        %  problems 
    end

    % ft = F;
    ft = y0t;
    % convert concentration
    AB_average = AB_average*vScale/v;
    F_average = F_average*vScale/v;

    % Reduce array size
    tt = [t(1:9) t(10:10:end)]*(10^xFactor);    % Time in days
    ft = [ft(1:9) ft(10:10:end)];
%     AB_average = [AB_average(1:9) AB_average(10:10:end)];
%     F_average = [F_average(1:9) F_average(10:10:end)];
    n10 = round((n10 - 9)/10);
    n90 = round((n90 - 9)/10);
    % n90 = find(tt>t(n90),1);
    % n10 = find(tt>t(n10),1);

    yLast = y(end);
    gLast = g(end);
    tLast = tt(end);
end % End AggKineticsV16PhaseExplorerCollection

function [testResult,testValue,n10,n90,x1,x2,AB_average,F_average] = TestAggValue(F,M,t,y0t)
    % Returns 1 if we have aggregation
    % there are two kinds of tests in this function
    % 1) compares avrage aggragation levels
    % 2) check for sustain aggragation by looking how many zero aggregation we have

    % for average aggregation levels
    % passvalratio = 75;      % Ratio higer than passvalratio indicates aggregation -> test fail
    % tPrecent = 0.1;         % precent of time from the begining or the end I want to check

    % The function always returns the average of the AB monomers at the end
    % of the test

    % For sustained non-zero aggregation level
    passvalratio = 0;        % Precent of zero aggregation  (consider less strict value like, sligthly larger than 0
    tPrecent = 0.20;         % precent of time from the begining or the end I want to check
    tmax = t(end);           % max time value

    % choosing a time 
    % Using the size of the time step as an indication to something happening
    lt = length(t);
    timesteps = (t(2:end) - t(1:end-1));
    mintimesteps = min(nonzeros(timesteps));
    timesteps = 3*timesteps/mintimesteps;
    timesteps = log(timesteps);
    firstpos = find(timesteps < (timesteps(1)/3));
    if isempty(firstpos)
        % startPos is the estimated aggregation onset time position
        startPos = 10;
        dn10 = 0;
        dn90 = 0;
    else
        startPos = firstpos(1);
        % the end location of initial aggregation levels
        dn10 = - round(0.2*startPos)*(firstpos(1)>20);           
        % when using aggregation level average
%         dn90 = round((lt*0.8-firstpos(1))*(firstpos(1)<lt*0.8) - 20*(firstpos(1)>(lt-20)));
        % when using aggregation sustainability do not use dn90
        % replace dn90 with zero
        dn90 = 0;
    end


    % I'm doing the avarage over 10% of the time
    % Force none-zero n10 and n90
    n10 = max([find(t>tPrecent*tmax,1),startPos+dn10,10]);
    n90 = max(find(t>(1-tPrecent)*tmax,1),startPos+dn90);
    n90 = min(n90,lt-11);

    % location for start aggregation sustainability test
    ns = round(2*n10);
    ns = ns*(ns<n90) + n90*(ns>n90);

    % avrage of acommulation over time at the begining and the end
    x1 = sum(F(2:n10).*(t(2:n10)-t(1:(n10-1))))/(t(n10)-t(1));
    % find the number no aggregation elements
    zf = length(find(F(ns:end)==0));
    zm = length(find(M(ns:end)==0));
    % x2 = sum(F(end-n90+1:end).*(t(end-n90+1:end)-t(end-n90:end-1)))/(t(end)-t(end-n90));
    x2 = sum(F(n90+1:end).*(t(n90+1:end)-t(n90:end-1)))/(t(end)-t(n90));
    % raltive acommulation 
    
    % for average test
    % testValue = x2/x1;
    % testResult = (testValue > passvalratio);
    
    % For sustained aggregation
    AB_average = sum(y0t(n90+1:end).*(t(n90+1:end)-t(n90:end-1)))/(t(end)-t(n90));
    F_average = x2;
    testValue = min([zf,zm])*100/(lt-ns+1);
    testResult = (testValue <= passvalratio);
end % end TestAggValue

function x = test_labs(DataArray)
    x=3;
end

function TransPos = GetTransionPos2(TestVal)
    % Get test transion position and direction
    % TransPosSign = +1 when going from 1 to 0
    % TransPosSign = -1 when going from 0 to 1
    % TransPosSign = 0 if all Test values are the same
    % [1,1,1,1] ->  -4
    % [0,0,0,0] ->  4
    % this function is used by ExplorPhaseSpace.m

    if length(TestVal)>4
        disp('TestVal array to long')
    end
    if max(TestVal)>1 || min(TestVal)<0
        disp('invalid values in the TestVal array')
    end
    % check transition position
    TransPosSign = ((TestVal(1)==0)-(TestVal(1)==1))*(sum(TestVal)>=0 && sum(TestVal)<=4);
    if (sum(TestVal)== 0) || (sum(TestVal)== 4)
        TransPos = 4;
    else
        TransPos = 1;
        while TestVal(TransPos)==TestVal(TransPos+1) && TransPos<3
            TransPos = TransPos + 1;
        end
    end

    TransPos = TransPos*TransPosSign;
end

function f = IsIn(x,a)

    % This function is x is in a
    % a would typically be an array

    f = 0;

    for i=1:length(a)
        if x==a(i)
            f = 1;
            break
        end
    end

end
