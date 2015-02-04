function f = IsIn(x,a)

% This function is x is in a
% a would typically be an array

f = 0;

for i=1:length(a)
    if x==a(i)
        f = 1;
        break
    end
end
