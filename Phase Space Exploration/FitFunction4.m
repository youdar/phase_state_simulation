function [nFit,fitParameter,fittedFun,regularSol,ExplorDir,r,baseAngle]=FitFunction4(i,x,y,xFactor,yFactor)

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


