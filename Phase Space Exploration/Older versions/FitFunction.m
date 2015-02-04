function [nFit,fitParameter,fittedFun,regularSol]=FitFunction(i,nFit,fitParameter,x,y,xFactor,yFactor)

% Find the best function to predict the next point
% do not use more that 6 points

yy = y/yFactor;
xx = x/xFactor;
testvalue = 0.97;

while ((i-nFit)>0 && fitParameter>=testvalue) && nFit<6
    % Start with linear fit
    [TempFun,gof]=fit(yy((i-nFit):i),xx((i-nFit):i),'poly1');
    fitParameter = abs(GetFitParameter(gof));
    % check if the line is close to horizontal, if yes use horizontal solution
    if sum(yy(i)==yy((i-nFit):i))~=(nFit+1)
        regularSol = 1;     
    else                    % we have a stright horizontal line  
        regularSol = 0;
    end
    
        
    %     if fitParameter<testvalue
%         % Try 2nd order polynomial fit
%         [TempFun,gof]=fit(yy((i-nFit):i),xx((i-nFit):i),'exp1');
%         fitParameter = abs(GetFitParameter(gof));
%     end
%     if fitParameter<testvalue
%         % Try 2nd order polynomial fit
%         [TempFun,gof]=fit(yy((i-nFit):i),xx((i-nFit):i),'poly2');
%         fitParameter = abs(GetFitParameter(gof));
%     end
%     if fitParameter<testvalue
%         % Try 'exp2' fit
%         [TempFun,gof]=fit(yy((i-nFit):i),xx((i-nFit):i),'exp2');
%         fitParameter = abs(GetFitParameter(gof));
%     end
    if fitParameter>=testvalue           % update the function and nFit
        nFit = nFit + 1;            % Add 1 to nFit if the linear fit is good
        fittedFun = TempFun;
    end
end

disp('Function fitted to data')
fittedFun
