function OutArray = InitParameterUpdate(initParameter,basePnum,NewY,adjPnum,NewX)

% This function adjusts initParameter array
% 1: klpf       elongation (monomer attachment) for free aggregates
% 2: klmf       monomer dettachment, for free aggregates 
% 3: klpm       elongation (monomer attachment) for membrane bounnd aggregates
% 4: klmm       monomer dettachment, for membrane bounnd aggregates
% 5: kap        attachemnt rate, of aggregates to te membrane
% 6: kam        dettachment rate, of aggregates to te membrane
% 7: kn         free aggregates nucleation rate
% 8: cnmRatio   cnmRatio = (cnm equivalent)/cn = cnm*p0/cn
% 9: kbf        free aggregates brakage rate
% 10: kbm       mambrane aggregates brakage rate
% 11: y00       initial A-Beta concentration
% 12: p00       initial Prion concentration
% 13: rs        Sequestration rate
% 14: rfm       Out-flow rate
% 15: n         size of maximum aggregates length
% 16: m         number of iterations
% 17: FailCounter

% Get the previous values
OutArray = initParameter;


% update values
flag1 = (basePnum == 11 || basePnum == 12);

% When keeping klmf and klmm equal
flag2 = (adjPnum == 2 || adjPnum == 4);
% When keeping kbf and  kbm equal
flag3 = (adjPnum == 9 || adjPnum == 10);
% When keeping klpf and klpm equal
flag4 = (adjPnum == 1 || adjPnum == 3);
% When keeping rs and rfm equal
flag5 = (adjPnum == 13 || adjPnum == 14);

if  flag1 && flag2 
    OutArray(basePnum) = NewY;
    OutArray(2) = NewX;
    OutArray(4) = NewX;
elseif flag1 && flag3
    OutArray(basePnum) = NewY;
    OutArray(9) = NewX;
    OutArray(10) = NewX;
elseif flag1 && flag4
    OutArray(basePnum) = NewY;
    OutArray(1) = NewX;
    OutArray(3) = NewX;
% elseif flag1 && flag5
%     OutArray(basePnum) = NewY;
%     OutArray(13) = NewX;
%     OutArray(14) = NewX;
else
    OutArray(basePnum) = NewY;
    OutArray(adjPnum) = NewX;
end