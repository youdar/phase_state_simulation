function [x1,y1,x2,y2,angleFactor] = SetTempXY(inputMatrix,TransPos,adjPnum,basePnum,x1,y1,x2,y2,OldTransPosSign)

% angleFactor indicates how the base angle should be adjusted, direction
% and units of angle resulotion at the particular run

switch TransPos         % done with this point - set values of y0,x0
    case -1             % [1,0,0,0] 
        x1 = inputMatrix(1,adjPnum);    % Set current value      
        y1 = inputMatrix(1,basePnum);   % Set current value  
        x2 = inputMatrix(2,adjPnum);    % Set current value      
        y2 = inputMatrix(2,basePnum);   % Set current value
        angleFactor = -1;
    case 1              % [0,1,1,1]
        x1 = inputMatrix(1,adjPnum);    % Set current value
        y1 = inputMatrix(1,basePnum);   % Set current value
        x2 = inputMatrix(1,adjPnum);    % Set current value
        y2 = inputMatrix(1,basePnum);   % Set current value
        angleFactor = -1;
    case -2             % [1,1,0,0]
        x1 = inputMatrix(2,adjPnum);    % Set current value      
        y1 = inputMatrix(2,basePnum);   % Set current value
        x2 = inputMatrix(3,adjPnum);    % Set current value      
        y2 = inputMatrix(3,basePnum);   % Set current value
        angleFactor = 0;
    case 2              % [0,0,1,1]
        x1 = inputMatrix(2,adjPnum);    % Set current value
        y1 = inputMatrix(2,basePnum);   % Set current value
        x2 = inputMatrix(2,adjPnum);    % Set current value
        y2 = inputMatrix(2,basePnum);   % Set current value
        angleFactor = 0;
    case -3             % [1,1,1,0]
        x1 = inputMatrix(3,adjPnum);    % Set current value
        y1 = inputMatrix(3,basePnum);   % Set current value
        x2 = inputMatrix(4,adjPnum);    % Set current value
        y2 = inputMatrix(4,basePnum);   % Set current value
        angleFactor = 1;
    case 3              % [0,0,0,1]
        x1 = inputMatrix(3,adjPnum);    % Set current value
        y1 = inputMatrix(3,basePnum);   % Set current value
        x2 = inputMatrix(4,adjPnum);    % Set current value
        y2 = inputMatrix(4,basePnum);   % Set current value
        angleFactor = 1;
    case -4             % [1,1,1,1]
        if OldTransPosSign >0
            % x1 stay the same
            % y1 stay the same
            x2 = inputMatrix(1,adjPnum);    % Set current value
            y2 = inputMatrix(1,basePnum);   % Set current value
            angleFactor = -2;
        elseif OldTransPosSign <0
            x1 = inputMatrix(4,adjPnum);    % Set current value
            y1 = inputMatrix(4,basePnum);   % Set current value
            % x2 stay the same
            % y2 stay the same
            angleFactor = 2;
        else
            disp('OldTransPosSign = 0   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            x1 = 0;
            y1 = 0;
            x2 = 0;
            y2 = 0;
            angleFactor = 0;
        end   
    case 4              % [0,0,0,0]
        if OldTransPosSign >0
            x1 = inputMatrix(4,adjPnum);    % Set current value
            y1 = inputMatrix(4,basePnum);   % Set current value
            % x2 stay the same
            % y2 stay the same 
            angleFactor = 2;
        elseif OldTransPosSign <0
            % x1 stay the same
            % y1 stay the same
            x2 = inputMatrix(1,adjPnum);    % Set current value
            y2 = inputMatrix(1,basePnum);   % Set current value
            angleFactor = -2;
        else
            disp('OldTransPosSign = 0   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            x1 = 0;
            y1 = 0;
            x2 = 0;
            y2 = 0;
            angleFactor = 0;
        end
end 

end
