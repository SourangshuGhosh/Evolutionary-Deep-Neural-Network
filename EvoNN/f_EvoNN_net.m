function out=f_EvoNN_net(data,Setslog)
%Setslog = importdata(Setslog);

nonodes = Setslog.nonodes;
noinnodes = Setslog.noinnodes;
nooutnodes = Setslog.nooutnodes;
f_in_index = Setslog.in_index;
f_out_index = Setslog.out_index;
xmin = Setslog.Data_min(1,f_in_index);
xmax = Setslog.Data_max(1,f_in_index);
ymin = Setslog.Data_min(1,f_out_index);
ymax = Setslog.Data_max(1,f_out_index);

if Setslog.subsets == 1
    i = Setslog.dataset.pareto.select;
    Net = Setslog.dataset.pareto.P(i);
else
    s = Setslog.subsets + 1;
    i = Setslog.dataset(s).pareto.select;
    Net = Setslog.dataset(s).pareto.P(i);
end

w = Net.w;
W = Net.W;

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
    for i = 1:nonodes,
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

