function plotparallel()
%A=importdata('cRVEAopt.mat'); % when not required xlsx remove %
[A,label]= xlsread('time2007obj1.xlsx',4);
roundn = @(x,n) round(x .* 10.^n)./10.^n; 
% B=A.FunctionValue;
 B=A(:,1:4);
 m=size(B);
n=m(2);
m=m(1);
for i=1:n
    if B(1,i)< 0
        M(:,i)=-B(:,i);
    else
        M(:,i)=B(:,i);
    end
 end
%label={'X1','X2','X3','X4','X5','X6','X7','X8','X9','X10','X11','X12','Y1','Y2','Y8'},
%label={'X1','X2','X3','X4','X5','X6','X7','X8','X9','X10','X11','X12'},
% label={'Y1','Y2','Y3','Y4'},
 label={'CO/CO2 (Y1)','O2 (Y2)','Productivity (Y3)','RAFT (Y4)'};
%label={'NG (X1)','PC (X2)','Coke Rate (Y1)','Carbon Rate (Y2)','Fuel Cost (Y3)'};
%%Figure out the lower bounds and ranges of the data...
Rgs = [max(M,[],1); min(M,[],1)];
nrm = Rgs(1,:)-Rgs(2,:);
shft = Rgs(2,:);
%%...so you can normalize them to [0 1]
Mnrm = (M - ones(m,1)*shft)./ (ones(m,1)*nrm);
%%Now plot the data, one line per row
figure(1);
plot(1:n, Mnrm,'Linestyle','-','Linewidth',2);
axis off; % switch off the axes
%%Because we're going to draw our own
for i=1:n
% A line for every column
line([i i],[0 1],'LineStyle','--','Color','k')
% Show the lower and upper end of the respective y-axis
text(i,0,num2str(roundn(Rgs(2,i),3)),'Fontsize',15);
text(i,1,num2str(roundn(Rgs(1,i),3)),'Fontsize',15);
text(i,1.05,label(i),'Color','r','Fontsize',20);
for j=0.25:0.25:0.75
    text(i,j,num2str(roundn(j*nrm(i)+shft(i),3)),'Fontsize',15);
end
end
%%Make the figure canvas white
set(gcf,'Color','w')
end
