function OutArray = InitUpdateVarArray2(initParameter,y,x,xFactor,yFactor,basePnum,adjPnum,ExplorDir,Extra_n,xDist,labNumber)
% This function is updating the values sent to each lab running the
% AggKineticsV10PhaseExplorer.m
% this function is used by ExplorPhaseSpace.m

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

OutArray = InitParameterUpdate(initParameter,basePnum,(y+dy),adjPnum,(x+dx));
OutArray(15) = round(OutArray(15)*Extra_n);        % Adjust n (aggregate length) if needed