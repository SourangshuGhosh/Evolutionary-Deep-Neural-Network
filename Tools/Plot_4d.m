clear all;
%Copy Pareto file to this directory
%Add legends, axes and title manually
%Save manually
[A,label]= xlsread('time2007obj1.xlsx',5);
%data = importdata('cRVEAopt.mat'); 
objectives = [1,2,3,4]; %Only 3 objectives
mrkSize = 200; %defines markersize
label={'CO/CO2 (Y1)','O2 (Y2)','Productivity (Y3)','RAFT (Y4)'};
%label={'Tuyere Velocity (Y3)','Heat Loss (Y4)','Corrected Productivity(Y5)','Coke Rate (Y6)'};
%label={'Corrected Productivity (Y5)','Coke Rate (Y6)','Plate Cooling Heat Loss (Y7)', 'Carbon Rate (Y8)'};
%label={'Coke Rate (Y1)','Carbon Rate (Y2)','Fuel Cost (Y3)'};

%=======================================
%FVal = data.FunctionValue(:,objectives);
B=A;
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
FVal = M(:,objectives);
%FVal = abs(FVal); %Use if you want to use absolute values of FVal
close all
colormap jet; % change this to change colour scheme
scatter3(FVal(:,1),FVal(:,2),FVal(:,3),mrkSize,FVal(:,4))
%scatter3(FVal(:,1),FVal(:,2),FVal(:,3), 'SizeData', 100, 'MarkerFaceColor','r')
xlabel(label(1),'FontSize',16, 'Rotation', 12, 'FontWeight', 'bold');
ylabel(label(2),'FontSize',16, 'Rotation', 340, 'FontWeight', 'bold');
zlabel(label(3),'FontSize',16, 'FontWeight', 'bold');
set(gca,'FontSize', 14, 'Box', 'on', 'LineWidth', 2);
colorbar