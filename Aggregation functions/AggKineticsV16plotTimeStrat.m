function [tOnset,i] = AggKineticsV16plotTimeStrat(t)
    % find the aggregation onset time
    
    % calculate dt
    timesteps = (t(2:end) - t(1:end-1));      
    % find the smallest non-zero dt
    mintimesteps = min(nonzeros(timesteps));       
    % normalized all dt using the smallest one. multiply by 3 to make sure
    % the log is larger than 1
    timesteps = 3*timesteps/mintimesteps;           
    timesteps = log(timesteps);
    
    firstpos = find(timesteps < (timesteps(1)/3),1);
    % This might give a better estimate for onset time, in some cases
%     i = find(timesteps < (timesteps(1)/2),1);
    if isempty(firstpos)
        disp('Problem with the aggregation onset time')
        startPos = 10;
        dn10 = 0;
    else
        startPos = firstpos;
        % Adjust the location of the 
        dn10 = round(0.2*startPos)*(startPos>20);           
    end
    i = startPos - dn10;
    tOnset = t(i);
end