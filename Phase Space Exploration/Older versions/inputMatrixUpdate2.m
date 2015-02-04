function inputMatrix = inputMatrixUpdate2(initParameter,y,dy,x,resFactor,basePnum,adjPnum,ExplorDir,Extra_n,xDist,nFit,xFactor,yFactor,fittedFun,regularSol)

% update and prepare the inputMatrix matrix for the 4 labs

if nFit>1
    newY = y+dy;
    yy = newY/yFactor;  
    dyy = dy/yFactor;                           % we evaluated the function with factorized arrays values
    if regularSol
        newX = feval(fittedFun,yy)*xFactor;     % get the x value
        dx = newX - x;
        baseAngle = atan(differentiate(fittedFun,yy))-pi/2; % the angle of the perpendicullar line at point y
    else
        dx = 0.1*x;
        newX = x + dx;
        baseAngle = -pi/2; % the angle of the perpendicullar line at point y
        if xDist==1
            xDist = 3;
        end
    end
    dxx = dx/xFactor;
    r = sqrt(dyy^2+dxx^2);
 
    inputMatrix(1,:) = UpdateVarArray3(initParameter,newY,newX,r,xFactor,yFactor,baseAngle,basePnum,adjPnum,Extra_n,xDist,resFactor,1);
    inputMatrix(2,:) = UpdateVarArray3(initParameter,newY,newX,r,xFactor,yFactor,baseAngle,basePnum,adjPnum,Extra_n,xDist,resFactor,2);
    inputMatrix(3,:) = UpdateVarArray3(initParameter,newY,newX,r,xFactor,yFactor,baseAngle,basePnum,adjPnum,Extra_n,xDist,resFactor,3);
    inputMatrix(4,:) = UpdateVarArray3(initParameter,newY,newX,r,xFactor,yFactor,baseAngle,basePnum,adjPnum,Extra_n,xDist,resFactor,4);                           
else        % for the case of 1st run with only one initial condition
    inputMatrix(1,:) = InitUpdateVarArray2(initParameter,y,x,xFactor,yFactor,basePnum,adjPnum,ExplorDir,Extra_n,xDist,1);
    inputMatrix(2,:) = InitUpdateVarArray2(initParameter,y,x,xFactor,yFactor,basePnum,adjPnum,ExplorDir,Extra_n,xDist,2);
    inputMatrix(3,:) = InitUpdateVarArray2(initParameter,y,x,xFactor,yFactor,basePnum,adjPnum,ExplorDir,Extra_n,xDist,3);
    inputMatrix(4,:) = InitUpdateVarArray2(initParameter,y,x,xFactor,yFactor,basePnum,adjPnum,ExplorDir,Extra_n,xDist,4);
end
  
% Displaying results on screen
disp(['+++++++++++++++   nFit = ' num2str(nFit) '  ++++++++++++++'])
for i=1:4
    disp(['y and x values: ' num2str(inputMatrix(i,basePnum)) ' , ' num2str(inputMatrix(i,adjPnum))])
end
if nFit>1
    disp(['r(' num2str(i) ')= ' num2str(r)])
    disp(['dx = ' num2str(dx)])
    disp(['dyy = ' num2str(dyy)])
    disp(['dxx = ' num2str(dxx)])
    disp(['yy = ' num2str(yy)])
    disp(['baseAngle = ' num2str(baseAngle*180/pi)])
end
    
