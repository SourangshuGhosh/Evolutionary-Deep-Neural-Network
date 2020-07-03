function r = sub(a,b)
%SUB Summary of this function goes here
%   Detailed explanation goes here

r = minus(a,b);
if r == zeros(size(a))
    r = inf*ones(size(a));
end

end