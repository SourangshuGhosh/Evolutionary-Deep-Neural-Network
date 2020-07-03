k=1;h=0;
for i= 1:3119
for j=i+1:3119
 if A(i,1)==A(j,1)&&A(i,2)==A(j,2)&&A(i,3)==A(j,3)&&A(i,4)==A(j,4)&&A(i,5)==A(j,5)&&A(i,6)==A(j,6)&&A(i,7)==A(j,7)&&A(i,8)==A(j,8)&&A(i,9)==A(j,9)&&A(i,10)==A(j,10)&&A(i,11)==A(j,11)&&A(i,12)==A(j,12)&&A(i,13)==A(j,13)&&A(i,14)==A(j,14)&&A(i,15)==A(j,15)&&A(i,16)==A(j,16)
     B(k)=j;
     k=k+1;
 end
end
[Y,PS] = removerows(A,'ind',B);
end
    