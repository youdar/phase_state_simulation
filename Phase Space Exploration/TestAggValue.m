function [testResult,testValue,n10,n90,x1,x2] = TestAggValue(F,M,t)
% Returns 1 is we have aggregation
% there are two kinds of tests in this function
% 1) compares avraage aggragation levels
% 2) check for sustain aggragation by looking how many zero aggregation we have

% for average aggregation levels
% passvalratio = 75;      % Ratio higer than passvalratio indicates aggregation -> test fail
% tPrecent = 0.1;         % precent of time from the begining or the end I
% want to check
% For sustained non-zero aggregation level
    passvalratio = 0;           % Precent of zero aggregation
    tPrecent = 0.20;         % precent of time from the begining or the end I want to check

    tmax = t(end);          % max time value



    % choosing a time 
    % Using the size of the time step as an indication to something happening
    lt = length(t);
    timesteps = (t(2:end) - t(1:end-1));
    mintimesteps = min(nonzeros(timesteps));
    timesteps = 3*timesteps/mintimesteps;
    timesteps = log(timesteps);
    firstpos = find(timesteps < (timesteps(1)/3));
    if isempty(firstpos)
        % startPos is the estimated aggregation onset time position
        startPos = 10;
        dn10 = 0;
        dn90 = 0;
    else
        startPos = firstpos(1);
        dn10 = - round(0.2*startPos)*(firstpos(1)>20);           % the end location of initial aggregation levels
        % when using aggregation level average
    %     dn90 = round((lt*0.8-firstpos(1))*(firstpos(1)<lt*0.8) - 20*(firstpos(1)>(lt-20)));
        % when using aggregation sustainability do not use dn90
        dn90 = 0;
    end


    % I'm doing the avarage over 10% of the time
    % Force none-zero n10 and n90
    n10 = max([find(t>tPrecent*tmax,1),startPos+dn10,10]);
    n90 = max(find(t>(1-tPrecent)*tmax,1),startPos+dn90);
    n90 = min(n90,lt-11);

    % location for start aggregation sustainability test
    ns = round(2*n10);
    ns = ns*(ns<n90) + n90*(ns>n90);

    % avrage of acommulation over time at the begining and the end
    x1 = sum(F(2:n10).*(t(2:n10)-t(1:(n10-1))))/(t(n10)-t(1));
    % fmfind the number no aggregation elements
    zf = length(find(F(ns:end)==0));
    zm = length(find(M(ns:end)==0));
    % x2 = sum(F(end-n90+1:end).*(t(end-n90+1:end)-t(end-n90:end-1)))/(t(end)-t(end-n90));
    x2 = sum(F(n90+1:end).*(t(n90+1:end)-t(n90:end-1)))/(t(end)-t(n90));
    % raltive acommulation 
    % for average test
    % testValue = x2/x1;
    % testResult = (testValue > passvalratio);
    % For sustained aggregation
    testValue = min([zf,zm])*100/(lt-ns+1);
    testResult = (testValue <= passvalratio);
end