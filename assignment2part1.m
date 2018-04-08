clc 
clear 
close all
 
load('LocVol.mat'); 
 
%% Numerical solution for interior points
Tint = T(2:length(T)-1); 
Kint = K(2:length(K)-1); 
for i=2:size(C,1)-1
   for j=2:size(C,2)-1
           dCdTint(i-1,j-1) = (C(i+1,j)-C(i-1,j))/(T(i+1)-T(i-1)); 
           dCdKint(i-1,j-1) = (C(i,j+1)-C(i,j-1))/(K(j+1)-K(j-1)); 
   end 
end 
 
deltaK = 1; 
for i=2:size(C,1)-1 
   for j=2:size(C,2)-1
           dCdK2int(i-1,j-1) = (C(i,j+1)-2*C(i,j)+C(i,j-1))/(deltaK)^2; 
   end 
end 
 
for i=1:size(dCdKint,1)
   for j=1:size(dCdKint,2)
      SigmaNumInt(i,j) = sqrt(2*(dCdTint(i,j)+(r-0)*Kint(j)*dCdKint(i,j)+0*C(i,j))/(Kint(j)*Kint(j)*dCdK2int(i,j))); 
   end 
end 
 
 
figure(1) 
[XnumInt,YnumInt]=meshgrid(Kint,Tint);
surf(XnumInt,YnumInt,SigmaNumInt);
xlabel("K")
ylabel("T")
zlabel("sigma")
title("Numerical solution for interior points")
%axis([60 200 0.5 1.5 0.1 0.3]);
axis tight
shading interp
colorbar

%% Numerical solution
%NEED TO ADD SQRT WHEN CALCULATING SIGMANUM. IT IS NOT ADDED BECAUSE A
%COUPLE OF POINTS ARE NEGATIVE, SO THE SQUARE ROOT GIVES AN ERROR. HOWEVER,
%IT WORKS USING ONLY INTERIOR POINTS
for i=1:size(C,1) 
   for j=1:size(C,2) 
       if i==1 && j==1 
           dCdT(i,j) = (C(i+1,j)-C(i,j))/(T(i+1)-T(i)); 
           dCdK(i,j) = (C(i,j+1)-C(i,j))/(K(j+1)-K(j)); 
       elseif i==1 && j<size(C,2) && j>1 
           dCdT(i,j) = (C(i+1,j)-C(i,j))/(T(i+1)-T(i)); 
           dCdK(i,j) = (C(i,j+1)-C(i,j-1))/(K(j+1)-K(j-1)); 
       elseif i==1 && j==size(C,2) 
           dCdT(i,j) = (C(i+1,j)-C(i,j))/(T(i+1)-T(i)); 
           dCdK(i,j) = (C(i,j)-C(i,j-1))/(K(j)-K(j-1)); 
       elseif i>1 && i<size(C,1) && j==1 
           dCdT(i,j) = (C(i+1,j)-C(i-1,j))/(T(i+1)-T(i-1)); 
           dCdK(i,j) = (C(i,j+1)-C(i,j))/(K(j+1)-K(j)); 
       elseif i>1 && i<size(C,1) && j==size(C,2) 
           dCdT(i,j) = (C(i+1,j)-C(i-1,j))/(T(i+1)-T(i-1)); 
           dCdK(i,j) = (C(i,j)-C(i,j-1))/(K(j)-K(j-1)); 
       elseif i==size(C,1) && j==1 
           dCdT(i,j) = (C(i,j)-C(i-1,j))/(T(i)-T(i-1)); 
           dCdK(i,j) = (C(i,j+1)-C(i,j))/(K(j+1)-K(j)); 
       elseif i==size(C,1) && j>1 && j<size(C,2) 
           dCdT(i,j) = (C(i,j)-C(i-1,j))/(T(i)-T(i-1)); 
           dCdK(i,j) = (C(i,j+1)-C(i,j-1))/(K(j+1)-K(j-1)); 
       elseif i==size(C,1) && j==size(C,2) 
           dCdT(i,j) = (C(i,j)-C(i-1,j))/(T(i)-T(i-1)); 
           dCdK(i,j) = (C(i,j)-C(i,j-1))/(K(j)-K(j-1)); 
       else 
           dCdT(i,j) = (C(i+1,j)-C(i-1,j))/(T(i+1)-T(i-1)); 
           dCdK(i,j) = (C(i,j+1)-C(i,j-1))/(K(j+1)-K(j-1)); 
       end 
   end 
end 
 
deltaK = 1; 
for i=1:size(C,1) 
   for j=1:size(C,2) 
       if j==1 
           dCdK2(i,j) = (2*C(i,j)-5*C(i,j+1)+4*C(i,j+2)-C(i,j+3))/((deltaK)^2); 
       elseif j==size(C,2) 
           dCdK2(i,j) = (2*C(i,j)-5*C(i,j-1)+4*C(i,j-2)-C(i,j-3))/((deltaK)^2); 
       else 
           dCdK2(i,j) = (C(i,j+1)-2*C(i,j)+C(i,j-1))/((deltaK)^2); 
       end 
   end 
end 
 
for i=1:size(dCdK,1) 
   for j=1:size(dCdK,2) 
      SigmaNum(i,j) = sqrt(abs(2*(dCdT(i,j)+(r-0)*K(j)*dCdK(i,j)+0*C(i,j))/(K(j)*K(j)*dCdK2(i,j)))); 
   end 
end
 
figure(2) 
[Xnum1,Ynum1]=meshgrid(K,T);
surf(Xnum1,Ynum1,SigmaNum);
xlabel("K")
ylabel("T")
zlabel("sigma")
title("Numerical solution")
%axis([60 200 0.5 1.5 0.1 0.3]);
axis tight
shading interp
colorbar


%% Numerical solution with implied volatilities
for i=1:size(C,1)
    for j=1:size(C,2)
        Timp = T(i);
        t = 0;
        
        myfun = @(sigmaImp,Cimp,Kimp,S0) normcdf(1/(sigmaImp*sqrt((Timp-t)))*(log(S0/Kimp)+(r-0+0.5*sigmaImp*sigmaImp)*(Timp-t)))*S0*exp(-0*(Timp-t))-normcdf(1/(sigmaImp*sqrt((Timp-t)))*(log(S0/Kimp)+(r-0-0.5*sigmaImp*sigmaImp)*(Timp-t)))*Kimp*exp(-r*(Timp-t))-Cimp;
        Kimp = K(j);
        Cimp = C(i,j);

        fun = @(sigmaImp) myfun(sigmaImp,Cimp,Kimp,S0);

        sigmaImp(i,j) = fzero(fun,0);
        d1(i,j) = 1/(sigmaImp(i,j)*sqrt(Timp))*(log(S0/Kimp)+(r-0+0.5*sigmaImp(i,j)*sigmaImp(i,j))*Timp);
        d2(i,j) = 1/(sigmaImp(i,j)*sqrt(Timp))*(log(S0/Kimp)+(r-0-0.5*sigmaImp(i,j)*sigmaImp(i,j))*Timp);
    end
end

for i=1:size(sigmaImp,1) 
   for j=1:size(sigmaImp,2) 
       if i==1 && j==1 
           dsdT(i,j) = (sigmaImp(i+1,j)-sigmaImp(i,j))/(T(i+1)-T(i)); 
           dsdK(i,j) = (sigmaImp(i,j+1)-sigmaImp(i,j))/(K(j+1)-K(j)); 
       elseif i==1 && j<size(sigmaImp,2) && j>1 
           dsdT(i,j) = (sigmaImp(i+1,j)-sigmaImp(i,j))/(T(i+1)-T(i)); 
           dsdK(i,j) = (sigmaImp(i,j+1)-sigmaImp(i,j-1))/(K(j+1)-K(j-1)); 
       elseif i==1 && j==size(sigmaImp,2) 
           dsdT(i,j) = (sigmaImp(i+1,j)-sigmaImp(i,j))/(T(i+1)-T(i)); 
           dsdK(i,j) = (sigmaImp(i,j)-sigmaImp(i,j-1))/(K(j)-K(j-1)); 
       elseif i>1 && i<size(sigmaImp,1) && j==1 
           dsdT(i,j) = (sigmaImp(i+1,j)-sigmaImp(i-1,j))/(T(i+1)-T(i-1)); 
           dsdK(i,j) = (sigmaImp(i,j+1)-sigmaImp(i,j))/(K(j+1)-K(j)); 
       elseif i>1 && i<size(sigmaImp,1) && j==size(sigmaImp,2) 
           dsdT(i,j) = (sigmaImp(i+1,j)-sigmaImp(i-1,j))/(T(i+1)-T(i-1)); 
           dsdK(i,j) = (sigmaImp(i,j)-sigmaImp(i,j-1))/(K(j)-K(j-1)); 
       elseif i==size(sigmaImp,1) && j==1 
           dsdT(i,j) = (sigmaImp(i,j)-sigmaImp(i-1,j))/(T(i)-T(i-1)); 
           dsdK(i,j) = (sigmaImp(i,j+1)-sigmaImp(i,j))/(K(j+1)-K(j)); 
       elseif i==size(sigmaImp,1) && j>1 && j<size(sigmaImp,2) 
           dsdT(i,j) = (sigmaImp(i,j)-sigmaImp(i-1,j))/(T(i)-T(i-1)); 
           dsdK(i,j) = (sigmaImp(i,j+1)-sigmaImp(i,j-1))/(K(j+1)-K(j-1)); 
       elseif i==size(sigmaImp,1) && j==size(sigmaImp,2) 
           dsdT(i,j) = (sigmaImp(i,j)-sigmaImp(i-1,j))/(T(i)-T(i-1)); 
           dsdK(i,j) = (sigmaImp(i,j)-sigmaImp(i,j-1))/(K(j)-K(j-1)); 
       else 
           dsdT(i,j) = (sigmaImp(i+1,j)-sigmaImp(i-1,j))/(T(i+1)-T(i-1)); 
           dsdK(i,j) = (sigmaImp(i,j+1)-sigmaImp(i,j-1))/(K(j+1)-K(j-1)); 
       end 
   end 
end 
 
deltaK = 1; 
for i=1:size(sigmaImp,1) 
   for j=1:size(sigmaImp,2) 
       if j==1 
           dsdK2(i,j) = (2*sigmaImp(i,j)-5*sigmaImp(i,j+1)+4*sigmaImp(i,j+2)-sigmaImp(i,j+3))/(deltaK)^2; 
       elseif j==size(C,2) 
           dsdK2(i,j) = (2*sigmaImp(i,j)-5*sigmaImp(i,j-1)+4*sigmaImp(i,j-2)-sigmaImp(i,j-3))/(deltaK)^2; 
       else 
           dsdK2(i,j) = (sigmaImp(i,j+1)-2*sigmaImp(i,j)+sigmaImp(i,j-1))/(deltaK)^2; 
       end 
   end 
end

for i=1:size(dsdK,1)
   for j=1:size(dsdK,2) 
      SigmaNumImp(i,j) = sqrt((sigmaImp(i,j)^2+2*sigmaImp(i,j)*T(i)*(dsdT(i,j)+r*K(j)*dsdK(i,j)))/(1+2*d1(i,j)*K(j)*sqrt(T(i))*dsdK(i,j)+K(j)*K(j)*T(i)*(d1(i,j)*d2(i,j)*dsdK(i,j)*dsdK(i,j)+sigmaImp(i,j)*dsdK2(i,j))));
   end 
end 
 
 
figure(3) 
[Xnum2,Ynum2]=meshgrid(K,T); 
surf(Xnum2,Ynum2,SigmaNumImp);
xlabel("K")
ylabel("T")
zlabel("sigma")
title("Numerical solution with implied volatilities")
axis([60 200 0.5 1.5 0.1 0.3]);
axis tight
shading interp
colorbar
 
%% Analitical solution 
figure(4)
[X,Y]=meshgrid(K,T); 
SigmaAnalyt = 0.15+0.15*(0.5+2.*Y).*((X./100-1.2).^2)./(((X.^2)./(100.^2))+1.44); 
surf(X,Y,SigmaAnalyt);
xlabel("K")
ylabel("T")
zlabel("sigma")
title("Analytical solution")
axis([60 200 0.5 1.5 0.1 0.3]);
axis tight
shading interp
colorbar
