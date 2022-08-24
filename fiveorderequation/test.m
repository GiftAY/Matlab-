clear all;
clc
A = [0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1; -243 -405 -270 -90 -15];
ts = 0;
tfinal = 1;
y0 = [0; 3; -9; -8; 0];
[y, t] = solution(@(t1,y)linearsystem(y, A), y0, 0, 5, 1e-5, 1e-5);
plot(t, y);
grid on
title('solution');
xlabel('x');
ylabel('t');
legend('t0', 't1', 't2', 't3', 't4');
figure
plot(t, y(:,1) - (-1/12)*exp(-3*t').*(129*t'.^4 + 16*t'.^3 - 54*t'.^2 - 36*t'))
grid on
title('difference');
xlabel('x');
ylabel('y-t0');