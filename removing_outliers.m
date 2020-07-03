clear;
[~,name,raw]=xlsread('jspl07out1.xls',1);
data1=xlsread('jsplobj072.xls',1);
datamirror=data1;
alloutliers=[];
for setn= 9:12
data=datamirror(:,setn);
means=mean(data);
strd=std(data);
withoutoutliers=[];
outliers=[];
x=2.5;
k=1;
l=1;
for i=1:length(data)
    if (means-x*strd) < data(i) && data(i) < (means+x*strd)
        withoutoutliers(k,:)=datamirror(i,:);
        k=k+1;
    end
    if means-x*strd > data(i)||data(i)> means+x*strd
        outliers(l,:)= datamirror(i,:);
        m(l)=i;
        l=l+1;
    end
end
datamirror=withoutoutliers;
alloutliers=[alloutliers;outliers];
end

%xlswrite('jspl08obj-no-outliers.xls',withoutoutliers);