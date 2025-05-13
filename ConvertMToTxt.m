load('frj712_5395_MoCoParam_7T_Series28.mat');
formatSpec = '%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f\n';
fileID = fopen('rp_adata.txt','w');
fprintf(fileID,formatSpec,R');
fclose(fileID);

%plot data
figure(1)
plot(R(1:end,1))

figure(2)
plot(R(1:end,2))

figure(3)
plot(R(1:end,3))

figure(4)
plot(R(1:end,4))

figure(5)
plot(R(1:end,5))

figure(6)
plot(R(1:end,6))