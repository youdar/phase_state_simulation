function OutArray = UpdateVarArray5(initParameter,y,x,r,baseAngle,xFactor,yFactor,basePnum,adjPnum,ExplorDir,Extra_n,resFactor,labNumber)

% This function is updating the values sent to each lab running the
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

OutArray = InitParameterUpdate(initParameter,basePnum,newY,adjPnum,newX);
OutArray(15) = round(OutArray(15)*Extra_n);        % Adjust n (aggregate length) if needed

end