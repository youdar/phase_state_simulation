function OutArray = UpdateVarArray3(initParameter,newy,newx,r,xFactor,yFactor,baseAngle,basePnum,adjPnum,Extra_n,xDist,resFactor,labNumber)
% This function is updating the values sent to each lab running the
% AggKineticsV10PhaseExplorer.m
% this function is used by ExplorPhaseSpace2.m

d = resFactor*0.01*r;

fArray0 = [-7.5,-2.5,2.5,7.5]*d;            % Initial search
fArray1 = [-1.5,-0.5,0.5,1.5]*d;            % Fine shearch, center region
fArray2 = [3.5,4.5,5.5,6.5]*d;              % Fine shearch, side regions
fArray3 = [-0.5,0,0,0.5]*d;                 % Horizontal line
switch abs(xDist)
    case 0     % First Cycle
        dyy = fArray0(labNumber)*sin(baseAngle);
        dxx = fArray0(labNumber)*cos(baseAngle);
    case 1
        dyy = fArray1(labNumber)*sin(baseAngle);
        dxx = fArray1(labNumber)*cos(baseAngle);
    case 2
        dyy = xDist/2*fArray2(labNumber)*sin(baseAngle);
        dxx = xDist/2*fArray2(labNumber)*cos(baseAngle);
    case 3
        dyy = fArray3(labNumber)*sin(baseAngle);
        dxx = fArray3(labNumber)*cos(baseAngle);
end

dy = dyy*yFactor;
dx = dxx*xFactor;

OutArray = InitParameterUpdate(initParameter,basePnum,(newy+dy),adjPnum,(newx+dx));
OutArray(15) = round(OutArray(15)*Extra_n);        % Adjust n (aggregate length) if needed
% OutArray(16) = round(OutArray(16)*ExtraIteration); % Adjust m (aggregate length) if needed