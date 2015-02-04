function [fittedFun,fitParameter] = FitFunctionDataPlot(x,y)

fitParameter = [0 0 0 0 0 0 0];

[TempFun1,gof1]=fit(x,y,'poly1');
fitParameter(1) = abs(GetFitParameter(gof1));
if fitParameter(1) <0.99
    [TempFun2,gof2]=fit(x,y,'poly2');
    fitParameter(2) = abs(GetFitParameter(gof2));
    [TempFun3,gof3]=fit(x,y,'poly3');
    fitParameter(3) = abs(GetFitParameter(gof3));
%     [TempFun4,gof4]=fit(x,y,'weibull');
%     fitParameter(4) = abs(GetFitParameter(gof4));
%     [TempFun5,gof5]=fit(x,y,'fourier4');
%     fitParameter(5) = abs(GetFitParameter(gof5));
    [TempFun6,gof6]=fit(x,y,'rat11');
    fitParameter(6) = abs(GetFitParameter(gof6));
    [TempFun7,gof7]=fit(x,y,'rat22');
    fitParameter(7) = abs(GetFitParameter(gof7));
    [TempFun8,gof8]=fit(x,y,'power2');
    fitParameter(8) = abs(GetFitParameter(gof8));
    [TempFun9,gof9]=fit(x,y,'exp1');
    fitParameter(9) = abs(GetFitParameter(gof9));
    [TempFun10,gof10]=fit(x,y,'exp2');
    fitParameter(10) = abs(GetFitParameter(gof10));
end

fitarray = fitParameter .* (fitParameter < 1);
% [val,ValPos] = max(fitarray);
ValPos = 8;
switch ValPos
    case 1
        fittedFun = TempFun1;
    case 2
        fittedFun = TempFun2;
    case 3
        fittedFun = TempFun3;
%     case 4
%         fittedFun = TempFun4;
%     case 5
%         fittedFun = TempFun5;
    case 6
        fittedFun = TempFun6;
    case 7
        fittedFun = TempFun7;
    case 8
        fittedFun = TempFun8;
    case 9
        fittedFun = TempFun9;
    case 10
        fittedFun = TempFun10;
end







