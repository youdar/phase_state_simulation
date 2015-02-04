function [x,y] = SetXY(inputMatrix,TransPos,adjPnum,basePnum,tempX,tempY,OldTransPosSign)

switch TransPos         % done with this point - set values of y0,x0
    case -1             % [1,0,0,0]
        x = inputMatrix(2,adjPnum);    % Set current value      
        y = inputMatrix(2,basePnum);   % Set current value   
    case 1              % [0,1,1,1]
        x = inputMatrix(1,adjPnum);    % Set current value
        y = inputMatrix(1,basePnum);   % Set current value
    case -2             % [1,1,0,0]
        x = inputMatrix(3,adjPnum);    % Set current value      
        y = inputMatrix(3,basePnum);   % Set current value
    case 2              % [0,0,1,1]
        x = inputMatrix(2,adjPnum);    % Set current value
        y = inputMatrix(2,basePnum);   % Set current value
    case -3             % [1,1,1,0]
        x = inputMatrix(4,adjPnum);    % Set current value
        y = inputMatrix(4,basePnum);   % Set current value
    case 3              % [0,0,0,1]
        x = inputMatrix(3,adjPnum);    % Set current value 
        y = inputMatrix(3,basePnum);   % Set current value
    case -4             % [1,1,1,1]
        x = tempX;
        y = tempY;
    case 4              % [0,0,0,0]
        if OldTransPosSign == 1
            x = inputMatrix(4,adjPnum);    % Set current value
            y = inputMatrix(4,basePnum);   % Set current value
        else
            x = inputMatrix(1,adjPnum);    % Set current value
            y = inputMatrix(1,basePnum);   % Set current value
        end
end 