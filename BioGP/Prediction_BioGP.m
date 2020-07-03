clear,clc
load sha.mat                 % load the new dataset in .mat format and not in cell form
load Setslog2-3-9_sha.mat     % load the Setslog file for the concerned objective
Data=num(80,:);
t=1:1:length(Data(:,1));
% y=Data(:,11);                         % Set the column number of the concerned objective
y=Data(:,Setslog.out_index);            % if new file is in same format column will be automatically taken
input_col=length(Setslog.in_index);
in=Data(:,1:input_col);
z=f_BioGP(in,Setslog);                 % z is the output
                                       % edit f_BioGP if the set is devided into subsets
plot(t,y,'*r',t,z,'-.k')
