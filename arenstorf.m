clear all;
clc;

y0 = [0.994, 0, 0, -2.00158510637908252240537862224]';
t0=0;
tfinal = 17.0652165601579625588917206249;
h = 0.0001;
tol = 1e-7;
tt = linspace(t0, tfinal, 100);


y0 = [1, 1, 1]';
t0 = 0;
tfinal = 20;
h = 0.001;
tol = 1e-9;
tt = linspace(t0, tfinal, 100);
[te1 ye1] = dormanprince(@lorenz, y0, t0, tfinal, h, tol);
figure
plot(ye1(1,:), ye1(2,:), 'b');
grid on

function dy = arenstorf(t,y)
mu = 0.012277471;
mu1 = 1 - mu;
D1 = ((y(1)+mu)^2 + y(2)^2)^(3/2);
D2 = ((y(1) - mu1)^2 + y(2)^2)^(3/2);

dy = [y(3); ...
      y(4); ...
      y(1) + 2*y(4) - mu1*(y(1)+mu)/D1 - mu*(y(1)-mu1)/D2; ...  
      y(2) - 2*y(3) - mu1*y(2)/D1 - mu*y(2)/D2;
      ];
end

function dy = lorenz(t,y)
sigma = 10;
beta = 8/3;
rho = 28;
dy = [sigma * (y(2) - y(1));
    y(1) * (rho - y(3)) - y(2);
    y(1) * y(2) - beta * y(3)];
end
