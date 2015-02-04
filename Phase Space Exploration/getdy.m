function dy = getdy(ExplorDir,y,dy1,yFactor,fittedFun,regularSol)

if regularSol
    dy = ExplorDir*min(abs(dy1),abs(2*y));      % Limit the size of the next step to twice the currnt value  
    dyy = dy/yFactor;
    yy = y/yFactor;

    if abs(atan(differentiate(fittedFun,yy+dyy)) - atan(differentiate(fittedFun,yy)))> pi/8
        dy = ExplorDir*0.001*y;                            % we evaluated the function with factorized arrays values 
    end
else
    dy = 0;
end