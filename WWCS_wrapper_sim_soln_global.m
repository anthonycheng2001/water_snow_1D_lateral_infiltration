%% Wrapper for WWCS_func

clear
clc
close all

%% Function Call

cold_temp_init = -1:1:-1;
%cold_temp_init = logspace(1,log10(273),10);
cold_temp = cold_temp_init;
N = length(cold_temp);
warm_temp = ones(1,N);
run_time = 50000*(1^2);

for i = 1:N

[phi,h,theta_w,theta_i,x,tdelt,T,Da,Ht,Pe,St] = WWCS_func_sym(warm_temp(i),cold_temp(i),run_time);

% timestamps
timestamp1 = round(run_time/1000);
timestamp2 = round(run_time/100);
timestamp3 = round(run_time/10);
timestamp4 = round(run_time);

% snapshots
snap1 = (timestamp1)*tdelt;
snap2 = (timestamp2)*tdelt;
snap3 = (timestamp3)*tdelt;
snap4 = (timestamp4)*tdelt;

% times [seconds]
time1 = snap1*T/60/60/24/365;
time2 = snap2*T/60/60/24/365;
time3 = snap3*T/60/60/24/365;
time4 = snap4*T/60/60/24/365;

% porosity
figure
plot(x, phi(1,:),'linewidth',2)
hold on
plot(x, phi(timestamp1,:),'linewidth',2)
plot(x,phi(timestamp2,:),'linewidth',2)
plot(x,phi(timestamp3,:),'linewidth',2)
plot(x, phi(timestamp4,:),'linewidth',2)
title(['Porosity for ' num2str(warm_temp(i)) '$^\circ$C Source, ' ...
    num2str(cold_temp(i)) '$^\circ$C Sink'],'FontSize',14,'Interpreter', ...
    'latex')
xlabel('Horizontal Position','FontSize',14,'Interpreter','latex')
ylabel('Porosity $\phi$','FontSize',14,'Interpreter','latex')
colororder([[4/5 1 1];[3/4 4/5 0.9];[2/3 3/4 0.8];[1/2 2/3 0.7];[0 1/2 0.6]])
legend('Initial',[num2str(time1) ' yrs after Initial'],[num2str(time2) ' yrs after Initial'],[num2str(time3) ' yrs after Initial'],[num2str(time4) ' yrs after Initial'],'fontsize',10, 'interpreter','latex','location','southeast')
hold off

% filename = sprintf('phi_w%0.1f_c%0.1f.mat',warm_temp(i),cold_temp(i));
% save(['/\jumbo\ice\infiltration\supraglcial_lake_snowplug_jess\porosity_data\' filename])
% filename = sprintf('phi_w%0.1f_c%0.1f.fig',warm_temp(i),cold_temp(i));
% savefig(filename)
% filename = sprintf('phi_w%0.1f_c%0.1f.png',warm_temp(i),cold_temp(i));
% saveas(gcf,filename)
% close all

filename = sprintf('phi_w%0.1f_c%0.1f.mat',warm_temp(i),cold_temp(i));
save(filename)
% filename = sprintf('phi_w%0.1f_c%0.1f.fig',warm_temp(i),cold_temp(i));
% savefig(filename)
filename = sprintf('phi_w%0.1f_c%0.1f.png',warm_temp(i),cold_temp(i));
saveas(gcf,filename)
close all

% water temperature
figure
plot(x, theta_w(1,:),'linewidth',2)
hold on
plot(x, theta_w(timestamp1,:),'linewidth',2)
plot(x,theta_w(timestamp2,:),'linewidth',2)
plot(x,theta_w(timestamp3,:),'linewidth',2)
plot(x, theta_w(timestamp4,:),'linewidth',2)
title(['Water Temperature for ' num2str(warm_temp(i)) '$^\circ$C Source, ' ...
    num2str(cold_temp(i)) '$^\circ$C Sink'],'FontSize',14,'Interpreter', ...
    'latex')
xlabel('Horizontal Position','FontSize',14,'Interpreter','latex')
ylabel('Water Temperature $\theta_w$','FontSize',14,'Interpreter','latex')
colororder([[4/5 1 1];[3/4 4/5 0.9];[2/3 3/4 0.8];[1/2 2/3 0.7];[0 1/2 0.6]])
legend('Initial',[num2str(time1) ' mins after Initial'],[num2str(time2) ' mins after Initial'],[num2str(time3) ' mins after Initial'],[num2str(time4) ' yrs after Initial'],'fontsize',10,'interpreter','latex')
hold off

filename = sprintf('theta_w_w%0.1f_c%0.1f.mat',warm_temp(i),cold_temp(i));
save(filename)
% filename = sprintf('theta_w_w%0.1f_c%0.1f.fig',warm_temp(i),cold_temp(i));
% savefig(filename)
filename = sprintf('theta_w_w%0.1f_c%0.1f.png',warm_temp(i),cold_temp(i));
saveas(gcf,filename)
close all

% ice temperature
figure
plot(x, theta_i(1,:),'linewidth',2)
hold on
plot(x, theta_i(timestamp1,:),'linewidth',2)
plot(x,theta_i(timestamp2,:),'linewidth',2)
plot(x,theta_i(timestamp3,:),'linewidth',2)
plot(x, theta_i(timestamp4,:),'linewidth',2)
title(['Ice Temperature for ' num2str(warm_temp(i)) '$^\circ$C Source, ' ...
    num2str(cold_temp(i)) '$^\circ$C Sink'],'FontSize',14,'Interpreter', ...
    'latex')
xlabel('Horizontal Position','FontSize',14,'Interpreter','latex')
ylabel('Snow Temperature $\theta_i$','FontSize',14,'Interpreter','latex')
colororder([[4/5 1 1];[3/4 4/5 0.9];[2/3 3/4 0.8];[1/2 2/3 0.7];[0 1/2 0.6]])
legend('Initial',[num2str(time1) ' mins after Initial'],[num2str(time2) ' mins after Initial'],[num2str(time3) ' mins after Initial'],[num2str(time4) ' yrs after Initial'],'fontsize',10,'interpreter','latex','location','southeast')
hold off

filename = sprintf('theta_i_w%0.1f_c%0.1f.mat',warm_temp(i),cold_temp(i));
save(filename)
% filename = sprintf('theta_i_w%0.1f_c%0.1f.fig',warm_temp(i),cold_temp(i));
% savefig(filename)
filename = sprintf('theta_i_w%0.1f_c%0.1f.png',warm_temp(i),cold_temp(i));
saveas(gcf,filename)
close all

% depth
figure
plot(x, h(1,:),'linewidth',2)
hold on
plot(x, h(timestamp1,:),'linewidth',2)
plot(x,h(timestamp2,:),'linewidth',2)
plot(x,h(timestamp3,:),'linewidth',2)
plot(x, h(timestamp4,:),'linewidth',2)
title(['Water Depth for ' num2str(warm_temp(i)) '$^\circ$C Source, ' ...
    num2str(cold_temp(i)) '$^\circ$C Sink'],'FontSize',14,'Interpreter', ...
    'latex')
xlabel('Horizontal Position','FontSize',14,'Interpreter','latex')
ylabel('Depth $h$','FontSize',14,'Interpreter','latex')
colororder([[4/5 1 1];[3/4 4/5 0.9];[2/3 3/4 0.8];[1/2 2/3 0.7];[0 1/2 0.6]])
legend('Initial',[num2str(time1) ' yrs after Initial'],[num2str(time2) ' yrs after Initial'],[num2str(time3) ' yrs after Initial'],[num2str(time4) ' yrs after Initial'],'fontsize',10,'interpreter','latex','location','southeast')
hold off

filename = sprintf('h_w%0.1f_c%0.1f.mat',warm_temp(i),cold_temp(i));
save(filename)
% filename = sprintf('h_w%0.1f_c%0.1f.fig',warm_temp(i),cold_temp(i));
% savefig(filename)
filename = sprintf('h_w%0.1f_c%0.1f.png',warm_temp(i),cold_temp(i));
saveas(gcf,filename)
close all

end

%% 0 source sink

% instantiate
t = linspace(0,time4,run_time+1);
phi_0 = 0.5;
tol = 1e-5;
W = warm_temp;
h_anal = tol*ones(size(h));
h_anal(:,1) = 1;
h_anal(:,end) = h_anal(:,end-1);

% solving for h
for j = 2:length(t)
    % h_anal(j,2:end-1) = h_anal(j-1,2:end-1)+Da.*phi_0.^2.*...
    %     tdelt./2./(x(2)-x(1)).^2.*(h_anal(j-1,3:end).*(h_anal(j-1,3:end)...
    %     -h_anal(j-1,2:end-1))-h_anal(j-1,1:end-2).*(h_anal(j-1,2:end-1)-...
    %     h_anal(j-1,1:end-2)));
    h_anal(j,2:end-1) = h_anal(j-1,2:end-1)+Da.*phi_0.^2.*tdelt./(x(2)-x(1)).^2.*...
        ((h_anal(j-1,3:end)-h_anal(j-1,1:end-2)).^2/4+h_anal(j-1,2:end-1).*...
        (h_anal(j-1,3:end)-2.*h_anal(j-1,2:end-1)+h_anal(j-1,1:end-2)));

    % flux on cold edge
    h_anal(j,1) = 1;
    h_anal(j,end) = h_anal(j,end-1);

    % crop
    h_anal(h_anal > 1) = 1;
    h_anal(h_anal < tol) = tol;

    disp('solving h analytically')
    disp(run_time-j)
end

figure
plot(x, h_anal(1,:),'linewidth',2)
hold on
plot(x, h_anal(timestamp1,:),'linewidth',2)
plot(x,h_anal(timestamp2,:),'linewidth',2)
plot(x,h_anal(timestamp3,:),'linewidth',2)
plot(x, h_anal(timestamp4,:),'linewidth',2)
title('Water Depth for $0^\circ$C source, $0^\circ$C Sink','FontSize',14,'Interpreter', ...
    'latex')
xlabel('Horizontal Position','FontSize',14,'Interpreter','latex')
ylabel('Depth $h$','FontSize',14,'Interpreter','latex')
colororder([[4/5 1 1];[3/4 4/5 0.9];[2/3 3/4 0.8];[1/2 2/3 0.7];[0 1/2 0.6]])
legend('Initial',[num2str(time1) ' yrs after Initial'],[num2str(time2) ' yrs after Initial'],[num2str(time3) ' yrs after Initial'],[num2str(time4) ' yrs after Initial'],'fontsize',10,'interpreter','latex','location','southeast')
hold off

filename = sprintf('h_w0_c0.mat');
save(filename)
filename = sprintf('h_w0_c0.png');
saveas(gcf,filename)
close all


%% check simulation for similarity

t = linspace(0,time4,run_time+1);
eta = zeros(length(t),length(x));
eta_var = eta;

PHI_s = 0.5;
d_p = 0.0004;
k_0 = PHI_s^2*d_p^2/180;
rho = 900;
g = 9.81;
mu = 1e-3;
phi_0 = 0.5;
gamma = k_0*rho*g/mu*phi_0^2;
gamma_var = k_0*rho*g./mu*phi.^2;
H = 1;
X = 260;

for j = 1:length(t)
    eta(j,:) = x./sqrt(gamma*H*T/X^2*t(j));
    eta_var(j,:) = x./sqrt(gamma_var(j,:)*H*T/X^2*t(j));
end

figure
plot(eta(2:end,1:5),h_anal(2:end,1:5),'linewidth',2)
hold on
title('Checking Similarity Assumption for $0^\circ$C','FontSize',14,'Interpreter','latex')
xlabel('$\eta = \frac{x}{\sqrt{\gamma\: t}}$','FontSize',12,'Interpreter','latex')
ylabel('$h(\eta)$','FontSize',12,'Interpreter','latex')
hold off

filename = sprintf('check_similarity_w0_c0.mat');
save(filename)
filename = sprintf('check_similarity_w0_c0.png');
saveas(gcf,filename)
close all

figure
plot(eta_var(2:end,1:5),h(2:end,1:5),'linewidth',2)
hold on
title('Checking Similarity Assumption for $\pm1^\circ$C','FontSize',25,'Interpreter','latex')
xlabel('$\eta = \frac{x}{\phi\sqrt{Da\: t}}$','FontSize',20,'Interpreter','latex')
ylabel('$h(\eta)$','FontSize',20,'Interpreter','latex')
hold off

filename = sprintf('check_similarity_w1_c-1.mat');
save(filename)
filename = sprintf('check_similarity_w1_c-1.png');
saveas(gcf,filename)
close all


%% shooting to solve f(eta)

% eta_N = 3;
% g_guess = 1;
% f = @(xf,yf)[-eta_N*yf(2); yf(3); -eta_N*(1-xf)/yf(2)/2*yf(3)-yf(3)^2/yf(2)];
% p1 = fzero(@(z)odetest1(g_guess,f,z),-1);
% p2 = fzero(@(z)odetest2(z,f),g_guess);
% yi = [p2;1e-3;p1];
% sol = ode15s(@(t,x)f(t,x),[1e-3,1],yi);
% f_end = sol.y(2,end);
% g_end = sol.y(1,end);
% g_flux_check = sol.y(1,1)+2*sol.y(3,end);

eta_end = 3;
f = @(xf,yf)[yf(2); -xf*yf(2)/2/yf(1)-yf(2)^2/yf(1)];
p1 = fzero(@(z)odetest1(z,f,eta_end),-1);
yi = [1;p1];
sol = ode15s(@(t,x)f(t,x),[0,eta_end],yi);
eta_N_index = find(sol.y(1,:)<1e-3,1,'first');
eta_N = sol.x(eta_N_index);

% 1, -1
p1_var = fzero(@(z)odetest1(z,f,eta_end),-1);
yi_var = [1;p1_var];
sol_var = ode15s(@(t,x)f(t,x),[0,eta_end],yi_var);
eta_N_index_var = find(sol_var.y(1,:)<1e-3,1,'first');
eta_N_var = sol_var.x(eta_N_index_var);

%% compare similarity approx to simulation similarity

figure
plot(sol.x,sol.y(1,:),'LineWidth',2)
hold on
plot(eta(2:end,1:5),h_anal(2:end,1:5),'k-','linewidth',2)
title('Global Constraint Shooting Method','FontSize',25,'Interpreter','latex')
xlabel('$\eta$','FontSize',20,'Interpreter','latex')
ylabel('$f(\eta)$','FontSize',20,'Interpreter','latex')
legend('Shooting Method','Simulation Data','location','northeast')
hold off

filename = sprintf('shoot_vs_simulation_global_w0_c0.mat');
save(filename)
filename = sprintf('shoot_vs_simulation_global_w0_c0.png');
saveas(gcf,filename)
close all

% 1,-1
figure
plot(sol_var.x,sol_var.y(1,:),'LineWidth',2)
hold on
plot(eta_var(2:end,1:5),h(2:end,1:5),'k-','linewidth',2)
title('Global Constraint Shooting Method','FontSize',25,'Interpreter','latex')
xlabel('$\eta$','FontSize',20,'Interpreter','latex')
ylabel('$f(\eta)$','FontSize',20,'Interpreter','latex')
legend('Shooting Method','Simulation Data','location','northeast')
hold off

filename = sprintf('shoot_vs_simulation_global_w1_c-1.mat');
save(filename)
filename = sprintf('shoot_vs_simulation_global_w1_c-1.png');
saveas(gcf,filename)
close all


%% reconstruct h

sol = ode15s(@(t,x)f(t,x),[0,eta_N],yi);

x_init = sol.x*sqrt(gamma*t(2));
x_1 = sol.x*sqrt(gamma*H*T/X^2*t(timestamp1));
x_2 = sol.x*sqrt(gamma*H*T/X^2*t(timestamp2));
x_3 = sol.x*sqrt(gamma*H*T/X^2*t(timestamp3));
x_4 = sol.x*sqrt(gamma*H*T/X^2*t(timestamp4));

% 1,-1
sol_var = ode15s(@(t,x)f(t,x),linspace(0,eta_N_var,2),yi_var);

x_init_var = sol_var.x.*sqrt(gamma_var(2,:).*t(2));
x_1_var = sol_var.x.*sqrt(gamma_var(timestamp1,:).*H.*T./X^2.*t(timestamp1));
x_2_var = sol_var.x.*sqrt(gamma_var(timestamp2,:).*H.*T./X^2.*t(timestamp2));
x_3_var = sol_var.x.*sqrt(gamma_var(timestamp3,:).*H.*T./X^2.*t(timestamp3));
x_4_var = sol_var.x.*sqrt(gamma_var(timestamp4,:).*H.*T./X^2.*t(timestamp4));

% x_init = sol.x*sqrt(t(2));
% x_1 = sol.x*sqrt(t(timestamp1));
% x_2 = sol.x*sqrt(t(timestamp2));
% x_3 = sol.x*sqrt(t(timestamp3));
% x_4 = sol.x*sqrt(t(timestamp4));

figure
plot(x_init,sol.y(1,:),'LineWidth',2)
hold on
plot(x_1,sol.y(1,:),'LineWidth',2)
plot(x_2,sol.y(1,:),'LineWidth',2)
plot(x_3,sol.y(1,:),'LineWidth',2)
plot(x_4,sol.y(1,:),'LineWidth',2)
plot(x, h_anal(1,:),'k--','linewidth',1)
plot(x, h_anal(timestamp1,:),'k--','linewidth',1)
plot(x,h_anal(timestamp2,:),'k--','linewidth',1)
plot(x,h_anal(timestamp3,:),'k--','linewidth',1)
plot(x, h_anal(timestamp4,:),'k--','linewidth',1)
title('Similarity Solution with Global Constraint','FontSize',25,'Interpreter', ...
    'latex')
xlabel('Horizontal Position','FontSize',20,'Interpreter','latex')
ylabel('Depth $h$','FontSize',20,'Interpreter','latex')
colororder([[4/5 1 1];[3/4 4/5 0.9];[2/3 3/4 0.8];[1/2 2/3 0.7];[0 1/2 0.6]])
legend('Initial',[num2str(time1) ' yrs after Initial'],[num2str(time2) ' yrs after Initial'],[num2str(time3) ' yrs after Initial'],[num2str(time4) ' yrs after Initial'],'Simulation','fontsize',10,'interpreter','latex','location','northeast')
hold off

filename = sprintf('h_similarity_global_w0_c0.mat');
save(filename)
filename = sprintf('h_similarity_global_w0_c0.png');
saveas(gcf,filename)
close all

%1,-1
figure
plot(x_init_var,sol_var.y(1,:),'LineWidth',2)
hold on
plot(x_1_var,sol_var.y(1,:),'LineWidth',2)
plot(x_2_var,sol_var.y(1,:),'LineWidth',2)
plot(x_3_var,sol_var.y(1,:),'LineWidth',2)
plot(x_4_var,sol_var.y(1,:),'LineWidth',2)
plot(x, h(1,:),'k--','linewidth',1)
plot(x, h(timestamp1,:),'k--','linewidth',1)
plot(x,h(timestamp2,:),'k--','linewidth',1)
plot(x,h(timestamp3,:),'k--','linewidth',1)
plot(x, h(timestamp4,:),'k--','linewidth',1)
title('Similarity Solution with Global Constraint for $\pm1^\circ$C','FontSize',25,'Interpreter', ...
    'latex')
xlabel('Horizontal Position','FontSize',20,'Interpreter','latex')
ylabel('Depth $h$','FontSize',20,'Interpreter','latex')
colororder([[4/5 1 1];[3/4 4/5 0.9];[2/3 3/4 0.8];[1/2 2/3 0.7];[0 1/2 0.6]])
legend('Initial',[num2str(time1) ' yrs after Initial'],[num2str(time2) ' yrs after Initial'],[num2str(time3) ' yrs after Initial'],[num2str(time4) ' yrs after Initial'],'Simulation','fontsize',10,'interpreter','latex','location','northeast')
hold off

filename = sprintf('h_similarity_global_w1_c-1.mat');
save(filename)
filename = sprintf('h_similarity_global_w1_c-1.png');
saveas(gcf,filename)
close all


%% functions

function c_out = odetest1(guess,fnctn,eta)
sol = ode15s(@(t,x)fnctn(t,x),[0,eta],[1;guess]);
%c_out = sol.y(1,1)+2*sol.y(3,end);
int = 0;
for j = 1:length(sol.x)-1
    int = int + 0.5*(sol.x(j+1)-sol.x(j))*(sol.y(1,j)+sol.y(1,j+1));
end
%c_out = sol.y(2,end)-1;
c_out = int+2*guess;
end

% function c_out = odetest2(gguess,fnctn)
% p1 = fzero(@(z)odetest1(gguess,fnctn,z),-1);
% sol = ode15s(@(t,x)fnctn(t,x),[1e-3,1],[gguess;1e-3;p1]);
% %c_out = sol.y(2,end)-1;
% c_out = sol.y(1,1)+2*sol.y(3,end);
% %c_out = sol.y(1,end);
% end

