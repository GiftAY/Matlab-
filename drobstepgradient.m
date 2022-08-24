clear all;
clc;
clear all;
clc;

%J = ||Ymes - X||^2 -> min(U)
%optimal conrol, bound condition, how will I use them?
% grJ_U = -2B'Ymes*tau + 2*tau*B'Ax + 2*BB'Uau^2 + 2*tau^2*B'C + 2*tau*B'x 

%constant parametres
lambda0 = 0.1;
eps = 0.05;
tau = 1;

A = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1];
B =  [1 0; 1 0; 1 0; 0 1; 0 0];
C = [1; 1; 1; -1; 1];
U0 = [-5*sin(2*pi/3); 10*cos(2*pi/3)];
Ymes = [179.43;    -404.38;    2.0407;    3.8453;       0]; %column
x = [179.11;    -404.23;    1.8492;    1.9926;    -4.2]; %column
xth = x;
kmax = 1000;
k = 1; i=2;
U01trace = U0(1);
U02trace = U0(2);
J_k = norm(Ymes - A*xth*tau - B*U0*tau - C*tau - xth)^2;
J_k1 = -10*J_k;

    while abs(J_k1 - J_k) > eps
        grJ_U = 2*B'*(-Ymes*tau + A*x*tau + B*U0*tau^2 + C*tau^2 + x*tau);
        U1 = U0;
        U0 = U0 - lambda0*grJ_U;
        J_k = norm(Ymes - A*xth*tau - B*U1*tau - C*tau - xth)^2;
        J_k1 = norm(Ymes - A*xth*tau - B*U0*tau - C*tau - xth)^2;
        %save the coordinates
        U01trace(i) = U0(1); 
        U02trace(i) = U0(2); 
        if( abs(J_k1 - J_k) < eps)
            l = 'accuracy';
            break
        end
        k = k+1;
        i = i+1;
    end    
U_opt = U0;
iterations_n = k;
%plot the soluion of he gradient descent 
[x, y] = meshgrid(-400:10:400, -100:10:100);
X = [x(:)'; y(:)'];
for i = 1:length(X)
    z(i) = (norm(Ymes - A*xth*tau - B*X(:,i)*tau - C*tau - xth))^2;
end
z=reshape(z,size(x));
grid,hold on
contour(x, y, z, '--');
plot(U01trace(:), U02trace(:), '-x');
%[C, h] = contour(x, y, z, '--');
%clabel(C,h);

