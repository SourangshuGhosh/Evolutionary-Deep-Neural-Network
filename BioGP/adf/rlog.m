function r = rlog(a)
%RLOG Summary of this function goes here
%   Detailed explanation goes here
if a > 0
    r = log(a);
else
    r = inf*ones(size(a));
end

end