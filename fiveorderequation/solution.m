function [y, t] = solution(f, y0, t, tfinal, h, tol)
facmax = 5; % the maximum of step increment
facmin = 0.2; %the minimum of step decrement
fac = 0.9; %fac is a garant factor
% Butcher's table of the coefficients
c2 = 1/5;  a21 = 1/5;
c3 = 3/10; a31 = 3/40; a32 = 9/40;
c4 = 4/5;  a41 = 44/45; a42 = -56/15; a43 = 32/9;
c5 = 8/9;  a51 = 19372/6561; a52 = -25360/2187; a53 = 64448/6561; a54 = -212/729;
c6 = 1;    a61 = 9017/3168; a62 = -355/33; a63 = 46732/5247; a64 = 49/176; a65 = -5103/18656;
c7 = 1;    a71 = 35/384; a72 = 0; a73 = 500/1113; a74 = 125/192; a75 = - 2187/6784; a76 = 11/84;
           b1 = 35/384; b2 = 0; b3 = 500/1113; b4 = 125/192; b5 = -2187/6784; b6 = 11/84; b7 =0;
           Z1 = 71/57600; Z2 =0; Z3 = -71/16695; Z4 = 71/1920; Z5 = -17253/339200; Z6 = 22/525; Z7 = -1/40;

y(1,:) = y0';
k = 2;
n = numel(y0);
while(t < tfinal)
    if (t+h) > tfinal
        h = tfinal - t;
    end
    yn = y(k-1,:)';
    tn = t(k-1);
    k1 = h*f(tn, yn);
    k2 = h*f(tn + c2*h, yn + a21*k1);
    k3 = h*f(tn + c3*h, yn + a31*k1 + a32*k2);
    k4 = h*f(tn + c4*h, yn + a41*k1 + a42*k2 + a43*k3);
    k5 = h*f(tn + c5*h, yn + a51*k1 + a52*k2 + a53*k3 + a54*k4);
    k6 = h*f(tn + c6*h, yn + a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5);
    k7 = h*f(tn + c7*h, yn + a71*k1 + a72*k2 + a73*k3 + a74*k4 + a75*k5 + a76*k6);
    y1 = yn + b1*k1 + b2*k2 + b3*k3 + b4*k4 + b5*k5 + b6*k6 + b7*k7;
    
    z1 = Z1*k1 + Z2*k2 + Z3*k3 + Z4*k4 + Z5*k5 + Z6*k6 + Z7*k7;
   
    err = abs(z1)./max([10^(-5)*ones(n, 1), abs(y0), abs(y1), 2*eps/tol*ones(n, 1)], [], 2);
    err = sqrt(norm(err)/n);
    
    h_new = h*min(facmax, max( facmin, fac*(tol/err)^0.2) );
    h = h_new;
    t(k) = t(k-1) + h;
    y(k,:) = y1;
    k = k+1;
end

