function [t, y] = domanprince(f, y0, t, tfinal, h, tol)  
                   
% Butcher's table of the coefficients
c2 = 1/5; a21 = 1/5;
c3 = 3/10; a31 = 3/40; a32 = 9/40;
c4 = 4/5;  a41 = 44/45; a42 = -56/15; a43 = 32/9;
c5 = 8/9; a51 = 19372/6561; a52 = -25360/2187; a53 = 64448/6561; a54 = -212/729;
c6 = 1; a61 = 9017/3168; a62 = -355/33; a63 = 46732/5247; a64 = 49/176; a65 = -5103/18656;
c7 = 1; a71 = 35/384; a72 = 0; a73 = 500/1113; a74 = 125/192; a75 = - 2187/6784; a76 = 11/84;
        b1 = 35/384; b2 = 0; b3 = 500/1113; b4 = 125/192; b5 = -2187/6784; b6 = 11/84; b7 =0;
        Z1 = 71/57600; Z2 =0; Z3 = -71/16695; Z4 = 71/1920; Z5 = -17253/339200; Z6 = 22/525; Z7 = -1/40;
% Dense output
        d1 = -1337/480; d2 = 1039/360; d3 = -1163/1152;
        d4 = 1054/9275; d5 = -4682/27825; d6 = 379/5565;
        d7 = 27/40; d8 = -9/5; d9 = 83/96; d10 = -3/250;
        d11 = 22/375; d12 = -37/600; d13 = -3/10; d14 = 29/30;
        d15 = -17/24;
% parametres
facmax = 5;
facmin = 0.1;
fac = 0.9;
h = 0.0001;
%iterations
k = 2;
tt = linspace(t,tfinal, 100);
n = numel(y0);
y(:,1) = y0
while(t < tfinal)
    tn = t(k-1);
    if(tn+h) > tfinal
        h = tfinal - tn;
    end
    yn = y(:,k-1);
    k1 = h * f(tn, yn);
    k2 = h * f(tn+h*c2, yn+k1*a21);
    k3 = h * f(tn+h*c3, yn+k1*a31 + k2*a32);
    k4 = h * f(tn+h*c4, yn + k1*a41 + k2*a42 + k3*a43);
    k5 = h * f(tn+h*c5, yn + k1*a51 + k2*a52 + k3*a53 + k4*a54);
    k6 = h * f(tn+h*c6, yn + k1*a61 + k2*a62 + k3*a63 + k4*a64 + k5*a65);
    k7 = h * f(tn+h*c7, yn + k1*a71+ k2*a72 + k3*a73 + k4*a74 + k5*a75 + k6*a76);

    Y = yn + b1*k1 + b2*k2 + b3*k3 + b4*k4  + b5*k5 + b6*k6 + b7*k7;
    z1 = Z1*k1 + Z2*k2 + Z3*k3 + Z4*k4 + Z5*k5 + Z6*k6 + Z7*k7;
    
    err = abs(z1)./max([10^(-5)*ones(n, 1), abs(y0), abs(Y), 2*eps/tol*ones(n, 1)], [], 2);
    err = sqrt(norm(err)/n);
    h_new = h*min( facmax, max( facmin, fac*(tol/err)^0.2) );
    h = h_new; 
    y(:,k) = Y;
    if err <= tol
         ppin = tt > tn & tt <= tn+h;
         if sum(ppin) > 0
         theta = (tt(ppin) - tn)/h;
         B = [theta.*(1+theta.*(d1 + theta.*(d2+theta.*d3))); ...
                0;...
                100*theta.^2.*(d4 + theta.*(d5 + theta.*d6))/3; ...
                -5*theta.^2.*(d7 + theta.*(d8 +theta.*d9))/2; ...
                18225*theta.^2.*(d10 +theta.*(d11 + theta.*d12))/848; ...
                -22*theta.^2.*(d13+theta.*(d14 + theta.*(d15)))/7 ];
            K = [k1, k2, k3, k4, k5, k6]';
            y(:,k-1) = yn' + B'*K;
         end
    end
    t(k) = t(k-1) + h; %3.285e-5 %5.223e-5 %6.731e-5
    k = k + 1;
end

end

        