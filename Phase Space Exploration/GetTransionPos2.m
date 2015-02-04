function TransPos = GetTransionPos2(TestVal)
% Get test transion position and direction
% TransPosSign = +1 when going from 1 to 0
% TransPosSign = -1 when going from 0 to 1
% TransPosSign = 0 if all Test values are the same
% [1,1,1,1] ->  -4
% [0,0,0,0] ->  4
% this function is used by ExplorPhaseSpace.m

if length(TestVal)>4
    disp('TestVal array to long')
end
if max(TestVal)>1 || min(TestVal)<0
    disp('invalid values in the TestVal array')
end
% check transition position
TransPosSign = ((TestVal(1)==0)-(TestVal(1)==1))*(sum(TestVal)>=0 && sum(TestVal)<=4);
if (sum(TestVal)== 0) || (sum(TestVal)== 4)
    TransPos = 4;
else
    TransPos = 1;
    while TestVal(TransPos)==TestVal(TransPos+1) && TransPos<3
        TransPos = TransPos + 1;
    end
end

TransPos = TransPos*TransPosSign;