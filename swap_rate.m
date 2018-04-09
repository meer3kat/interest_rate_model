%% Assignment 2 part 2 Calibration of a intereste rate model from STIBOR and swap rates
%
%
% In this part, we will calibrate the simply_compounded interest rates and swap 
% rates using the Vasicek model.
%% 1 calibrate paremeter 
% In this part, we calibrae the parameters $\theta$ = {k, $\phi$, 
% $\sigma$,r0}.by calling fminsearch and fminunc.

clear all;


x0 = [0.9572 0.4854 0.8003 0.1419];
%x0 = [1 1 1 1];
T1 = [1 7 30 60 90 180]; %in days
T1 = T1/360; % in years
Lmarket = [0.0300 0.0302 0.0307 0.0314 0.0320 0.0335];

T2 = [1 2 3 4 5 6 7 8 9 10];%in years swap maturities time
Smarket = [0.0351 0.0368 0.0375 0.0379 0.0381 0.0383 0.0384 0.0385 0.0385 0.0386];
swaptime = 0.25;
fun = @(x)cirmodel(x,T1,T2,Lmarket,Smarket,swaptime);
x1 = fminsearch(fun,x0); %from fminsearch
x2 = fminunc(fun, x0); %from fminunc

for i = 1:1:size(T1,2)
    tau(i) = T1(i) - 0;
    B(i) = (1 - exp(-x1(2)* tau(i)))/x1(2);
    A(i) = (x1(1) - (x1(3)*x1(3)/(2*x1(2)*x1(2)))) * (B(i) - tau(i)) - (x1(3)*x1(3)*B(i)*B(i)/(4*x1(2)));
    Z(i) = exp(A(i) - B(i)*x1(4));
    L(i) = (1 - Z(i))/(tau(i)*Z(i));
    tau2(i) = T1(i) - 0;
    B2(i) = (1 - exp(-x2(2)* tau(i)))/x2(2);
    A2(i) = (x2(1) - (x2(3)*x2(3)/(2*x2(2)*x2(2)))) * (B(i) - tau(i)) - (x2(3)*x2(3)*B(i)*B(i)/(4*x2(2)));
    Z2(i) = exp(A(i) - B(i)*x2(4));
    L2(i) = (1 - Z(i))/(tau(i)*Z(i));
    
    
end

for j = 1:size(T2,2)
    
    tau(j) = T2(j) - 0;
    B(j) = (1 - exp(-x1(2)* tau(j)))/x1(2);
    A(j) = (x1(1) - (x1(3)*x1(3)/(2*x1(2)*x1(2)))) * (B(j) - tau(j)) - (x1(3)*x1(3)*B(j)*B(j)/(4*x1(2)));
    Z(j) = exp(A(j) - B(j)*x1(4));
    n = 1;
    for k = swaptime:swaptime:T2(j)
        
        Bj(n) = (1 - exp(-x1(2)*k))/x1(2);
        Aj(n) = (x1(1) - (x1(3)*x1(3)/(2*x1(2)*x1(2)))) * (Bj(n) - k) - (x1(3)*x1(3)*Bj(n)*Bj(n)/(4*x1(2)));
        Zj(n) = exp(Aj(n) - Bj(n)*x1(4));
        n = n + 1;
    
    end
    
    S(j) = (1 - Z(j))/(swaptime * sum(Zj));
    
    tau2(j) = T2(j) - 0;
    B2(j) = (1 - exp(-x2(2)* tau(j)))/x2(2);
    A2(j) = (x2(1) - (x2(3)*x2(3)/(2*x2(2)*x2(2)))) * (B(j) - tau(j)) - (x2(3)*x2(3)*B(j)*B(j)/(4*x2(2)));
    Z2(j) = exp(A(j) - B(j)*x2(4));
    n = 1;
    for k = swaptime:swaptime:T2(j)
        
        Bj2(n) = (1 - exp(-x2(2)*k))/x2(2);
        Aj2(n) = (x2(1) - (x2(3)*x2(3)/(2*x2(2)*x2(2)))) * (Bj(n) - k) - (x2(3)*x2(3)*Bj(n)*Bj(n)/(4*x2(2)));
        Zj2(n) = exp(Aj(n) - Bj(n)*x2(4));
        n = n + 1;
    
    end
    
    S2(j) = (1 - Z(j))/(swaptime * sum(Zj));
    
end

my_fig = figure('position', [0, 0, 700, 500]);

subplot(2,1,1);
plot(T1*360,L*100);
hold on
plot(T1*360,L2*100,'g:');
plot(T1*360, Lmarket*100,'rd','MarkerSize',10);
xlim([0 T1(end)*360]);
ylim([2 5]);
xlabel('time(days)');
ylabel('yield(%)');
legend('interest rate(fminsearch)','interest rate(fminunc)','market data')
title('yield curve for the simply-compounded interest rates');
hold off

subplot(2,1,2); 
plot(T2,S*100);
hold on
plot(T2,S2*100,'g:');
plot(T2, Smarket*100,'rd','MarkerSize',10);
xlim([0 T2(end)]);
ylim([2 5]);
xlabel('time(year)');
ylabel('yield(%)');
legend('swap rate(fminsearch)','swap rate(fminunc)','market data')
title('yield curve for the swap rates');
hold off



%saveas(my_fig,'yield_plot2.png');


%% comment
% the different x0 will give different results when call fminsearch and
% fminunc. But when calculate back the simply-compouded interest rate and
% swap rates, they give the same result. This might be that the fminsearch
% and fminunc found local minimizer and the value are both close to 0.
% the convergency rate are fast, but fminsearch use Nelder_Mead simplex
% method that does not require the minimizing problem to be differentiable.
% the error function can be discontinuous. while the fminunc use the calssic newton method assuming that it is
% differentiable. fminsearch tend to give a more robust result while
% fminunc might get stucked at local minimum.

dbtype('cirmodel.m');


