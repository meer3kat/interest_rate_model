
x0 = [0.9572 0.4854 0.8003 0.1419];
x1 = fminsearch(@cirmodel,x0)
x2 = fminunc(@cirmodel, x0)

x1 =

    0.8625    0.4640    1.1628    0.0295


x2 =

    0.8630    0.2964    0.8933    0.0297