clear,clc
load Test_data.mat                     % load the new dataset in .mat and it must not be in cell form
load Setslog5-3_data_demo.mat     % load the Setslog file for the concerned objective
Data=Test_data(:,:);
t=1:1:length(Data(:,1));
% y=Data(:,11);                         % Set the column number of the concerned objective
y=Data(:,Setslog.out_index);     % or it will be taken from the dataset
%plot(t,y);
in=Data(:,Setslog.in_index);
z=f_EvoNN_net(in,Setslog);                    % z is the output from the neural network
% plot(t,z)
plot(t,y,'*r',t,z,'-.k')

