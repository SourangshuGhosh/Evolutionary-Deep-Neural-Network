function r = rsqrt(a)
%rsqrt Summary of this function goes here
%   Detailed explanation goes here
if a > 0
    r = sqrt(a);
else
    r = inf*ones(size(a));
end

end