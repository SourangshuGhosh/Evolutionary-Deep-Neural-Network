function out=evaluate_obj(data,Setslog)
%Setslog = importdata(Setslog);

nonodes = Setslog.nonodes;
noinnodes = Setslog.noinnodes;
nooutnodes = Setslog.nooutnodes;
xmin = Setslog.xmin;
xmax = Setslog.xmax;
ymin = Setslog.ymin;
ymax = Setslog.ymax;
w = Setslog.w;
W = Setslog.W;

x= data(:,1:noinnodes);

%scale the inputs 
Xmin = Setslog.Xmin;
Xmax = Setslog.Xmax;
xsc=[];
% xsc=Setslog.DataSet_sc;
for i=1:noinnodes
    xsc=[xsc Xmin+ (x(:,i)-xmin(i))/(xmax(i)-xmin(i))*(Xmax-Xmin)];
end
in = xsc;
    
    noexp = length(in(:,1));
    s = zeros(noexp,1); z = [];
    for i = 1:nonodes
        s(:) = w(i,1)*ones(noexp,1)+((w(i,2:noinnodes+1)*in')');
        z(:,i) = 1./(1+exp(-s(:)));
    end
    A = [ones(noexp,1) z];

    bber = A*W;

% rescale the outputs 
for i=1:nooutnodes
    bber(:,i)=ymin(i)+(ymax(i)-ymin(i))*(bber(:,i)-Xmin)/(Xmax-Xmin);
end

out=bber;

