function r = powsq(a)
%pow2 Summary of this function goes here
%   Detailed explanation goes here

r = power(a,2);
if r == zeros(size(a))
    r = inf*ones(size(a));
end

end