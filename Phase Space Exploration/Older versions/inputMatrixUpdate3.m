function inputMatrix = inputMatrixUpdate3(initParameter,y,dy,x,resFactor,basePnum,adjPnum,ExplorDir,Extra_n,xDist,nFit,xFactor,yFactor,fittedFun,regularSol)

% update and prepare the inputMatrix matrix for the 4 labs

inputMatrix(1,:) = UpdateVarArray4(initParameter,y,dy,x,xFactor,yFactor,basePnum,adjPnum,ExplorDir,Extra_n,xDist,nFit,fittedFun,resFactor,regularSol,1);
inputMatrix(2,:) = UpdateVarArray4(initParameter,y,dy,x,xFactor,yFactor,basePnum,adjPnum,ExplorDir,Extra_n,xDist,nFit,fittedFun,resFactor,regularSol,2);
inputMatrix(3,:) = UpdateVarArray4(initParameter,y,dy,x,xFactor,yFactor,basePnum,adjPnum,ExplorDir,Extra_n,xDist,nFit,fittedFun,resFactor,regularSol,3);
inputMatrix(4,:) = UpdateVarArray4(initParameter,y,dy,x,xFactor,yFactor,basePnum,adjPnum,ExplorDir,Extra_n,xDist,nFit,fittedFun,resFactor,regularSol,4);                           
  
% Displaying results on screen
disp(['+++++++++++++++   nFit = ' num2str(nFit) '  ++++++++++++++'])
for i=1:4
    disp(['y and x values: ' num2str(inputMatrix(i,basePnum)) ' , ' num2str(inputMatrix(i,adjPnum))])
end