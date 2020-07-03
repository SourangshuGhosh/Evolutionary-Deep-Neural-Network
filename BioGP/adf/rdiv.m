function r = rdiv(a,b)
%RDIV Summary of this function goes here
%   Detailed explanation goes here

r = rdivide(a,b);
if r == zeros(size(a))
    r = inf*ones(size(a));
end    

end