function OutArray = UpdateVarArray4(initParameter,y,dy,x,xFactor,yFactor,basePnum,adjPnum,ExplorDir,Extra_n,xDist,nFit,fittedFun,resFactor,regularSol,labNumber)
% This function is updating the values sent to each lab running the
% AggKineticsV10PhaseExplorer.m
% this function is used by ExplorPhaseSpace.m

if nFit>1 
    dyy = dy/yFactor;                                   % we evaluated the function with factorized arrays values
    if regularSol
        newX = feval(fittedFun,(y+dy)/yFactor)*xFactor;             % get the x value
        dx = newX - x;
        baseAngle = atan(differentiate(fittedFun,y/yFactor));  % the angle of the perpendicullar line at point y
        
    else        % Horizontal solution
        dx = 0.025*x;
        baseAngle = 0; % the angle of the perpendicullar line at point y
        if xDist==1
            xDist = 3;
        end
    end
    dxx = dx/xFactor;
    r = sqrt(dyy^2+dxx^2);
    fArray0 = [-7.5,-2.5,2.5,7.5]*resFactor;            % Initial search
    fArray1 = [-1.5,-0.5,0.5,1.5]*resFactor;            % Fine shearch, center region
    fArray2 = [3.5,4.5,5.5,6.5]*resFactor;              % Fine shearch, side regions
    fArray3 = [-0.5,0,0,0.5]*resFactor;                 % Horizontal line
    % Because we fitted x=f(y) the x and y axis are switched
    switch abs(xDist)
        case 0      % First Cycle
            dyy = r*cos(baseAngle+fArray0(labNumber));
            dxx = r*sin(baseAngle+fArray0(labNumber));
        case 1      % 2nd cycle, TransPos->{2,-2}
            dyy = r*cos(baseAngle+fArray1(labNumber));
            dxx = r*sin(baseAngle+fArray1(labNumber));
        case 2      % 2nd cycle, xDist=2,TransPos->{3,-3},xDist=-2,TransPos->{1,-1]}
            dyy = r*cos(baseAngle+xDist/2*fArray2(labNumber));  
            dxx = r*sin(baseAngle+xDist/2*fArray2(labNumber));  
        case 3      % 2nd cycle
            dyy = r*cos(baseAngle+fArray3(labNumber));
            dxx = r*sin(baseAngle+fArray3(labNumber));   
    end
    dy = ExplorDir*dyy*yFactor;
    dx = ExplorDir*dxx*xFactor;
else                % nFit = 1  -> Only one initial point
    r = ExplorDir*abs(0.025*y/yFactor);
    switch abs(xDist)
        case 0      % First Cycle
            baseAngle = 0;
            dy = r*sin(baseAngle+(labNumber-1)*pi/4);
            dx = r*cos(baseAngle+(labNumber-1)*pi/4);
        case 1      % Before 2nd cycle, TransPos = 1
            baseAngle = 0;
            dy = r*sin(baseAngle+labNumber*pi/20);
            dx = r*cos(baseAngle+labNumber*pi/20);
        case 2      % Before 2nd cycle, TransPos = 2
            baseAngle = pi/4;
            dy = r*sin(baseAngle+labNumber*pi/20);
            dx = r*cos(baseAngle+labNumber*pi/20);
        case 3      % Before 2nd cycle, TransPos = 3
            baseAngle = pi/2;
            dy = r*sin(baseAngle+labNumber*pi/20);
            dx = r*cos(baseAngle+labNumber*pi/20);
        case 4      % Before 2nd cycle, TransPos = 4
            baseAngle = 3*pi/4;
            dy = r*sin(baseAngle+labNumber*pi/20);
            dx = r*cos(baseAngle+labNumber*pi/20);   
    end
    dy = dy*yFactor;
    dx = dx*xFactor; 
end
newY = y+dy;
newX = x+dx;

OutArray = InitParameterUpdate(initParameter,basePnum,newY,adjPnum,newX);
OutArray(15) = round(OutArray(15)*Extra_n);        % Adjust n (aggregate length) if needed