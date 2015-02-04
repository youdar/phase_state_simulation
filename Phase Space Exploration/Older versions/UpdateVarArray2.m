function OutArray = UpdateVarArray2(initParameter,newy,newx,r,xFactor,yFactor,baseAngle,basePnum,adjPnum,Extra_n,xDist,resFactor,labNumber)
% This function is updating the values sent to each lab running the
% AggKineticsV10PhaseExplorer.m
% this function is used by ExplorPhaseSpace2.m

d = resFactor*0.01*r;
fArray1 = [-5.5,-0.5,0.5,5.5]*d;          % Search spread
fArray2 = [1,2,3,4]*d;              % Search spread

if xDist==0     % First Cycle
    dyy = fArray1(labNumber)*sin(baseAngle);
    dxx = fArray1(labNumber)*cos(baseAngle);
else            % xDist is -1 or 1
    dyy = xDist*fArray2(labNumber)*sin(baseAngle);
    dxx = xDist*fArray2(labNumber)*cos(baseAngle);
end

dy = dyy*yFactor;
dx = dxx*xFactor;

OutArray = InitParameterUpdate(initParameter,basePnum,(newy+dy),adjPnum,(newx+dx));
OutArray(15) = round(OutArray(15)*Extra_n);        % Adjust n (aggregate length) if needed
% OutArray(16) = round(OutArray(16)*ExtraIteration); % Adjust m (aggregate length) if needed