function y_out = cirmodel(x)
%here x is vector input, 
%x(1) = phi
%x(2) = k
%x(3) = sigma (volatility)
%x(4) = r0(present interest rate)

T1 = [1 7 30 60 90 180]; %in days
T1 = T1/360;
Lmarket = [0.0300 0.0302 0.0307 0.0314 0.0320 0.0335];%,market libor rate

%B = (1 - exp(-x(2)* tauti))/x(2);%need to get tauti

for i = 1:size(T1)
    tau(i) = T1(i) - 0;
    B(i) = (1 - exp(-x(2)* tau(i)))/x(2);
    A(i) = (x(1) - (x(3)*x(3)/(2*x(2)*x(2)))) * (B(i) - tau(i)) - (x(3)*x(3)*B(i)*B(i)/(4*x(2)));
    Z(i) = exp(A(i) - B(i)*x(4));
    L(i) = (1 - Z(i))/(tau(i)*Z(i));
    Ldiff2(i) = ((L(i) - Lmarket(i))/Lmarket(i))^2;
end

T2 = [1 2 3 4 5 6 7 8 9 10];%in years swap maturities time
Smarket = [0.0351 0.0368 0.0375 0.0379 0.0381 0.0383 0.0384 0.0385 0.0385 0.0386];

swaptime = 3/12;

for j = 1:size(T2)
    
    tau(j) = T2(j) - 0;
    B(j) = (1 - exp(-x(2)* tau(j)))/x(2);
    A(j) = (x(1) - (x(3)*x(3)/(2*x(2)*x(2)))) * (B(j) - tau(j)) - (x(3)*x(3)*B(j)*B(j)/(4*x(2)));
    Z(j) = exp(A(j) - B(j)*x(4));
    n = 1;
    for k = swaptime:swaptime:T2(j)
        
        Bj(n) = (1 - exp(-x(2)*k))/x(2);
        Aj(n) = (x(1) - (x(3)*x(3)/(2*x(2)*x(2)))) * (Bj(n) - k) - (x(3)*x(3)*Bj(n)*Bj(n)/(4*x(2)));
        Zj(n) = exp(Aj(n) - Bj(n)*x(4));
        n = n + 1;
    
    end
    
    S(j) = (1 - Z(j))/(swaptime * sum(Zj));
    Sdiff2(j) = ((S(j) - Smarket(j))/Smarket(j))^2;
end

y_out = sum(Ldiff2) + sum(Sdiff2);
end

      








