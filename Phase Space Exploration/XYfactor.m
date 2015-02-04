function [xFactor,yFactor] = XYfactor(x,y)

    xFactor = 1;
    yFactor = 1;
    % factorizing the x and y so that the fitting function will work better
    if x~=0
        xFactor = 10^floor(log10(x));      
    end
    if y~=0
        yFactor = 10^floor(log10(y));
    end
end