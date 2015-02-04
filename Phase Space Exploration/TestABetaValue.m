function [testResult,testValue,n10,n90,x1,x2] = TestABetaValue(y1,t)
% Testing the change in A-Beta concetration due to aggregate formation
% 
% this function test the upper limit of the transition. it does not explor
% whether there is a range of values for which the desired ration is kept

passvalratio = 60/100;       % the wanted test value

tmax = t(end);          % max time value
tPrecent = 0.05;         % precent of time from the begining or the end I want to check

% choosing a time 
% Using the size of the time step as an indication to something happening
lt = length(t);
timesteps = (t(2:end) - t(1:end-1));
mintimesteps = min(nonzeros(timesteps));
timesteps = 3*timesteps/mintimesteps;
timesteps = log(timesteps);
firstpos = find(timesteps < (timesteps(1)/3));
if isempty(firstpos)
    startPos = 10;
    dn90 = 0;
    dn10 = 0;
else
    startPos = firstpos(1);
    dn10 = - round(0.2*startPos)*(firstpos(1)>20);           % the end location of initial aggregation levels
    dn90 = round((lt*0.8-firstpos(1))*(firstpos(1)<lt*0.8) - 20*(firstpos(1)>(lt-20)));
end



% % I'm doing the avarage over 10% of the time
% n10 = find(t>tPrecent*tmax,1);
% n90 = find(t>(1-tPrecent)*tmax,1);

% % avrage of acommulation over time at the begining and the end
% x1 = sum(y1(2:n10).*(t(2:n10)-t(1:(n10-1))))/(t(n10)-t(1));
% x2 = sum(y1(end-n90+1:end).*(t(end-n90+1:end)-t(end-n90:end-1)))/(t(end)-t(end-n90));

% Force none-zero n10 and n90
n10 = max([find(t>tPrecent*tmax,1),startPos+dn10,10]);
n90 = max(find(t>(1-tPrecent)*tmax,1),startPos+dn90);
n90 = min(n90,lt-11);

% avrage of acommulation over time at the begining and the end
% x1 = sum(y1(2:n10).*(t(2:n10)-t(1:(n10-1))))/(t(n10)-t(1));
x2 = sum(y1(n90+1:end).*(t(n90+1:end)-t(n90:end-1)))/(t(end)-t(n90));
x1 = y1(1);
% x2 = sum(y1(n90:end))/(lt-n90+1);



testValue = x2/x1;
testResult = (testValue < passvalratio);

end