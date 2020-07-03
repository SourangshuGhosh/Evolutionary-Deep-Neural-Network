function [fval,W,InfoC] = NNevalnet(w)
global Setslog setno F_bad
global nonodes noinnodes nooutnodes

in = Setslog.dataset(setno).in;
out = Setslog.dataset(setno).out;
noexp = (length(out(:,1)));
Xmin = Setslog.Xmin;
Xmax = Setslog.Xmax;
ymin = Setslog.Data_min(1,Setslog.out_index);
ymax = Setslog.Data_max(1,Setslog.out_index);

s = zeros(noexp,1);
for i = 1:nonodes,
    s = w(i,1)*ones(noexp,1)+((w(i,2:noinnodes+1)*in')');
    z(:,i) = 1./(1+exp(-s));
end
A = [ones(noexp,1) z];
b = out;

%Solving the linear least square root problem
warning off
W = A\b;
bber = A*W;
% fval = sqrt(sum(((bber-b)/(max(b)-min(b))).^2)/noexp);
% add the AIC, AICc and BIC calculation here itself 

% rescale the outputs 
for i = 1:nooutnodes
    bber(:,i) = ymin(i)+(ymax(i)-ymin(i))*(bber(:,i)-Xmin)/(Xmax-Xmin);
    y(:,i) = ymin(i)+(ymax(i)-ymin(i))*(b(:,i)-Xmin)/(Xmax-Xmin);
end
fval = sqrt(sum(((bber-y)/(ymax-ymin)).^2)/noexp);

if isnan(fval) || isinf(fval)
    fval = F_bad+eps;
end
if ~isempty(find(bber > 1, 1)) || ~isempty(find(bber < 0, 1))
    fval = fval+0.05*(length(find(out < 0)) + length(find(out > 1)));
end

k=length(find(w(:,1:noinnodes+1)))+length(find(W));
n = noexp;
rss = sum(sum((bber-y).^2));

% InfoC.AIC = 2*k+n*log(rss/n);
% InfoC.AICc = InfoC.AIC+(2*k*(k+1)/(n-k-1));
% InfoC.BIC = k*log(k)+n*log(rss/n);
% InfoC.k = k;
% InfoC.n = n;
% InfoC.rss = rss;

InfoC = 2*k+n*log(rss/n);
InfoC = InfoC+(2*k*(k+1)/(n-k-1));

warning('on');
end