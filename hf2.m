clc
close all
clear all
format long
syms fi_ r_ a_ s_ l_ fi2_ A_ t_

%Data:
rho = 7800;
r = 0.08;
l = 0.06;
l_min = 0.01;
l_max = 0.2;
a = 0.02;
m = 1;
s = 600;

x_dir = r_*(1-cos(fi_))+a_;
y_dir = r_*sin(fi_);

len = (x_dir^2+y_dir^2)^0.5;
delta = len-l_;
F_fi = s_*delta;
U_fi = 0.5*s_*delta^2;
U_num = subs(subs(subs(subs(U_fi,l_,l),s_,s),a_,a),r_,r);
eq_eq0 = subs(subs(len,a_,a),r_,r)==l_;
eq_eq = subs(eq_eq0,l_,l);
eq_zero = solve(eq_eq);

dU = diff(U_num);
eq_eq_p = dU==0;
eq_points = solve(eq_eq_p)

% functions for plots
fi2 = linspace(-1*pi(),1*pi(),37);
U_0 = subs(U_num,fi_,fi2_);
fi_dot = ((4/(m*r^2))*(U_0-U_num)).^0.5;

% create wariables for potential plot
fi_draw = linspace(-2*pi(),2*pi(),361);
U_draw = subs(U_num,fi_,fi_draw);
x_draw = subs(subs(subs(s_*(x_dir-l),s_,s),a_,a),r_,r)*(-cos(fi_));
y_draw = subs(subs(s_*(y_dir),s_,s),r_,r)*(sin(fi_));

% create wariables for bifurcation plot
range_var = 200;
l_range = linspace(l_min,l_max,range_var);
l_range_2 = linspace(0.1285,l_max,range_var/2);
fi_bif = zeros(range_var,2);
fi_bif_2 = zeros(range_var/2,2);
for i=1:range_var
    eq_bif = subs(eq_eq0,l_,l_range(i));
    fi_bif(i,:) = solve(eq_bif);
end
for i=1:range_var/2
    eq_bif_2 = subs(eq_eq0,l_,l_range_2(i));
    fi_bif_2(i,:) = solve(eq_bif_2);
end

% task 4
l_new = subs(subs(len,a_,0.1),r_,r);
f_eq = sin(fi_)*(2*s/(m*r*l_new))*(l_new-l)*(0.1+r);
f_eq_dot = diff(f_eq)
f_eq_dotdot = diff(f_eq_dot)
f_eq_dotdotdot = diff(f_eq_dotdot)
omega_n_2 = eval(subs(f_eq_dot,fi_,0))
nu = (1/6)*eval(subs(f_eq_dotdotdot,fi_,0))
beta = (omega_n_2+3*nu*((A_)^2)/4)^0.5;
T_n = 2*pi()/beta;

x_0_t = A_*cos(beta*t_);
x_1_t = ((A_^3)/(32*beta^2))*(-cos(beta*t_)+cos(3*beta*t_));
fi_t = x_0_t+nu*x_1_t;
%% 
syms f_0_ omega_

omega_vals = linspace(0, 5000, 5000);
n_points = length(omega_vals);
f_0_all = [0, 0.5, 1];

% Preallocate: store all 3 roots per omega_plot
A_roots = NaN(3,n_points, 3);
for j = 1:3
    f_0 = f_0_all(j);
    for i = 1:n_points
        omega = omega_vals(i);
    
        % Coefficients of cubic: a*A^3 + b*A + c = 0
        a = 0.75 * nu;
        b = omega_n_2 - omega;
        c = -f_0 * omega_n_2;
    
        coeffs = [a 0 b c];
    
        % Solve cubic
        roots_i = roots(coeffs);
    
        % Store all 3 roots (real or complex)
        A_roots(j,i, :) = roots_i.';
    end
end
%% plot for task 2
fig = figure(1);
tiledlayout(2,1);
nexttile(1);
hold on
plot(linspace(eq_zero(1),eq_zero(1),2),linspace(0,max(U_draw),2),'--','Color',"#808080");
plot(linspace(eq_zero(2),eq_zero(2),2),linspace(0,max(U_draw),2),'--','Color',"#808080");
plot(linspace(2*pi()+eq_zero(1),2*pi()+eq_zero(1),2),linspace(0,max(U_draw),2),'--','Color',"#808080");
plot(linspace(-2*pi()+eq_zero(2),-2*pi()+eq_zero(2),2),linspace(0,max(U_draw),2),'--','Color',"#808080");
plot(fi_draw,U_draw)
hold off
grid on
grid minor
xlim([-2*pi(),2*pi()])
xticks(-2*pi():(pi/4):2*pi())
xlabel('$\varphi [rad]$','interpreter','latex')
ylabel('$U(\varphi) [J]$','interpreter','latex')

nexttile(2)
hold on
for i = 1:36
    color = rand(1, 3);
    plot(fi_draw,real(subs(subs(fi_dot,fi2_,fi2(i)),fi_,fi_draw)),'Color', color)
    plot(fi_draw,-real(subs(subs(fi_dot,fi2_,fi2(i)),fi_,fi_draw)), 'Color', color)
end
plot(linspace(eq_zero(1),eq_zero(1),2),linspace(-60,60,2),'--','Color',"#808080");
plot(linspace(eq_zero(2),eq_zero(2),2),linspace(-60,60,2),'--','Color',"#808080");
plot(linspace(2*pi()+eq_zero(1),2*pi()+eq_zero(1),2),linspace(-60,60,2),'--','Color',"#808080");
plot(linspace(-2*pi()+eq_zero(2),-2*pi()+eq_zero(2),2),linspace(-60,60,2),'--','Color',"#808080");
plot(linspace(-2*pi(),2*pi(),2),zeros(2),'Color',"#808080",'LineWidth',1);
%plot(zeros(2),linspace(-60,60,2),'Color',"#808080");

grid on
grid minor
xlim([-2*pi(),2*pi()])
xticks(-2*pi():(pi/4):2*pi())
xlabel('$\varphi [rad]$','interpreter','latex')
ylabel('$\dot{\varphi} [rad/s]$','interpreter','latex')
hold off
%f.Position = [700,300,800,1000];

%% plot for task 3
fig2 = figure(2);
hold on
plot(l_range,real(fi_bif(:,2)),"b")
plot(l_range_2,real(fi_bif_2(:,1)),"b")
plot(l_range,-real(fi_bif(:,2)),"b")
plot(l_range_2,-real(fi_bif_2(:,1)),"b")
plot(linspace(l_range_2(1),l_range_2(1),2),linspace(-fi_bif_2(1,1),fi_bif_2(1,1),2),'--','Color',"b");
plot(linspace(l_min,(a+2*r),2),linspace(-pi(),-pi(),2),'-','Color',"r");
plot(linspace(l_min,(a+2*r),2),linspace(pi(),pi(),2),'-','Color',"r");
plot(linspace(a,l_max,2),linspace(0,0,2),'-','Color',"r");
grid on
grid minor
yticks(-1.5*pi():(pi/4):1.5*pi())
ylabel('$\varphi [rad]$','interpreter','latex')
xlabel('$l [m]$','interpreter','latex')
hold off
%f.Position = [700,300,800,1000];
%% plot for task 4
% f_fi
fi_taylor = linspace(-2,2,100);
f_plot = subs(f_eq,fi_,fi_taylor);
taylor = subs(f_eq,fi_,0)+fi_*subs(f_eq_dot,fi_,0)+0.5*(fi_^2)*subs(f_eq_dotdot,fi_,0)+(1/6)*(fi_^3)*subs(f_eq_dotdotdot,fi_,0);
f_taylor = subs(taylor,fi_,fi_taylor);

A = linspace(0,pi/2,100);
T_n_plot = subs(T_n,A_,A);
T = linspace(0,100,1000);
fi_t_vec = subs(fi_t,A_,1);
fi_t_plot = double(subs(fi_t_vec,t_,T));

% Range of initial angles (in radians)
phi0_vals = linspace(0.001, pi()/2 - 0.001, 50);
periods = zeros(size(phi0_vals));

% Time span for each simulation
tspan = [0 10];

for i = 1:50
    phi0 = phi0_vals(i);
    x0 = [phi0, 0];  % zero initial velocity

    % Solve the ODE
    [t, x] = ode45(@phi_system, tspan, x0);

    phi = x(:,1);

    % Find zero crossings (from negative to positive)
    crossing_idx = find(phi(1:end-1) < 0 & phi(2:end) > 0);

    if length(crossing_idx) >= 2
        t_crossings = t(crossing_idx) - phi(crossing_idx) .* ...
            (t(crossing_idx+1) - t(crossing_idx)) ./ ...
            (phi(crossing_idx+1) - phi(crossing_idx));
        % One full period between two consecutive crossings
        periods(i) = t_crossings(2) - t_crossings(1);
    else
        periods(i) = NaN;  % Not enough oscillation detected
    end
end


fig3 = figure(3);
%tiledlayout(2,1);
%nexttile(1);
hold on
plot(fi_taylor,f_taylor,"b")
plot(fi_taylor,f_plot,"k")
hold off
xlabel('$\varphi [rad]$','interpreter','latex');
ylabel('$f(\varphi) [1/s^2]$','interpreter','latex');
legend({'Approximation','Original function','f_0 = 1'},'Location','northwest','FontSize',10)
ylim([-1000,1000])
grid on;
grid minor
%nexttile(2);
%plot(T,fi_t_plot)

% Amplitude
fig4 = figure(4);
hold on
plot(A,real(T_n_plot),"r")
plot(phi0_vals, periods, 'b-');
hold off
xlabel('$A [rad]$','interpreter','latex');
ylabel('$T [s]$','interpreter','latex');
legend({'Approximation','ode45','f_0 = 1'},'Location','northeast','FontSize',10)
xlim([0,pi()/10])
grid on;
grid minor

% Amplitude close
fig5 = figure(5);
hold on
plot(A,real(T_n_plot),"r")
plot(phi0_vals, periods, 'b-');
hold off
xlabel('$A [rad]$','interpreter','latex');
ylabel('$T [s]$','interpreter','latex');
legend({'Approximation','ode45','f_0 = 1'},'Location','northeast','FontSize',10)
grid on;
grid minor

%% 

fig6 = figure(6);
hold on
plot(omega_vals,real(A_roots(1,:,1)),"r")
plot(omega_vals(1082:end),abs(real(A_roots(2,1082:end,1))),"b")
plot(omega_vals(1082:end),abs(real(A_roots(3,1082:end,1))),"g")
plot(omega_vals,real(A_roots(1,:,2)),"r")

plot(omega_vals(2213:end),abs(real(A_roots(2,2213:end,2))),"b")
plot(omega_vals(2213:end),abs(real(A_roots(2,2213:end,3))),"b")
plot(omega_vals(1:1078),abs(real(A_roots(2,1:1078,3))),"b")


plot(omega_vals(2879:end),abs(real(A_roots(3,2879:end,2))),"g")
plot(omega_vals(2879:end),abs(real(A_roots(3,2879:end,3))),"g")
plot(omega_vals(1:1078),abs(real(A_roots(3,1:1078,3))),"g")

hold off
legend({'f_0 = 0','f_0 = 0.5','f_0 = 1'},'Location','northwest','FontSize',10)
xlabel('$\omega [1/s]$','interpreter','latex');
ylabel('$A [rad]$','interpreter','latex');
grid on;
grid minor

function dxdt = phi_system(t, x)
    % x(1) = phi, x(2) = phi_dot
    fi = x(1);
    fi_dot = x(2);
    
    r = 0.08;
    l = 0.06;
    m = 1;
    s = 600;
    a_=0.1;
    x_dir = r*(1-cos(fi))+a_;
    y_dir = r*sin(fi);
    
    l_new = (x_dir^2+y_dir^2)^0.5;

    % Compute angular acceleration
    fi_ddot = -sin(fi)*(2*s/(m*r*l_new))*(l_new-l)*(0.1+r);

    % Return derivatives
    dxdt = [fi_dot; fi_ddot];
end


