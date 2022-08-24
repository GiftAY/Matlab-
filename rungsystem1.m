clear all
grid, hold on

% psi_w = x(1) w = x(2) r = x(3) 

u_down = -6; psi_0 = -16.7;

%u_down = -4.332; psi_0 = -12.6; %u_down = -4.32; psi_0 = -12.89;

r_0 = 1; w_0 = 2 ;
u_tilda = psi_0/(2*r_0^2);
u_up = 100;



if (u_tilda > u_up)
    u = u_up;
elseif (u_tilda < u_down)
    u = u_down;
else
    u = u_tilda;
end


u0=u;

f=@(t,x)[-2*x(1)-4*(1-x(2)^2)*x(2)*x(3)^2; 2*x(2)+u; -x(3)];
Mesh = 10:5:100;

for k = 1:length(Mesh)
    N=Mesh(k); h=1/(N-1);
    y(:,1) = [psi_0, w_0, r_0];
    u_res = [u0];
        for n=1:(N-1)
            x(n)=(n-1)*h;
            k1=f(x(n),y(:,n));
            k2=f(x(n)+0.5*h,y(:,n)+0.5*h*k1);
            k3=f(x(n)+0.5*h,y(:,n)+0.5*h*k2);
            k4=f(x(n)+h,y(:,n)+h*k3);
            y(:,n+1) = y(:,n)+(h/6)*(k1+2*k2+2*k3+k4);
            
            u_tilda = y(1,n+1)/(2*y(3,n+1)^2);
            
            if (u_tilda > u_up)
                u = u_up;
            elseif (u_tilda < u_down)
                u = u_down;
            else 
                u = u_tilda;
            end
                u_res = [u_res, u];
        end
        
end
figure
grid, hold on
plot(0:h:1, u_res, '-');
ylabel('value of u')
xlabel('The coordinate of t')


figure
grid, hold on
plot(0:h:1,  y(1,:));






            
           


