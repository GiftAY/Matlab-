close all

clear
syms t psi(t) w(t) w(t) 

psi0 = sym('-12.402');
w0 = sym('2');
r0 = sym('1');
psi = psi0;
r = r0;
w = w0;

u_up = sym('7');
u_down = sym('-7');
u= sym('-6.201');

u_tilda = psi0/(2*r0^2);

iteration = [-2*psi-4*(1-w^2)*w*r^2; 2*w+u; -r];
nmax = 7;
for i=1:nmax
    iteration = [psi0; w0; r0] + int(iteration, 0, t);
    psi=iteration(1);
    w=iteration(2);
    r=iteration(3);
    
    u_tilda = psi/(2*r^2);
    a = u_tilda;
   
    iteration = [-2*psi-4*(1-w^2)*w*r^2; 2*w+u; -r];
end
figure
fplot(psi,[0 1]), grid, hold on
fplot(w,[0 1])
fplot(r,[0 1])
ylim([-15 20])


