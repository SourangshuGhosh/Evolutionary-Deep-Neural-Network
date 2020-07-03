function r = mul(a,b)
%MUL Summary of this function goes here
%   Detailed explanation goes here

r = times(a,b);
if r == zeros(size(a))
    r = inf*ones(size(a));
end

end
