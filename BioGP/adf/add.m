function r = add(a,b)
%ADD Summary of this function goes here
%   Detailed explanation goes here

r = plus(a,b);
if r == zeros(size(a))
    r = inf*ones(size(a));
end

end