%% Wrapper for WWCS_func

clear
clc
close all

%% Function Call

% solver setup as an array of sink temps
cold_temp_init = -5:1:-1;
cold_temp = flip(cold_temp_init);

% keep source at melting point for each sink temp run
N = length(cold_temp);
warm_temp = 0*ones(1,N);

% run time of ~1.8 years, visibly near steady state
run_time = 50000;

% solve for each round
for i = 1:N

    % function call
    [phi,h,theta_w,theta_i,x,tdelt,T,Da,Ht,Pe,St] = WWCS_func_rev6(warm_temp(i),cold_temp(i),run_time);

    % timestamps in run_time clicks
    timestamp1 = round(run_time/1000);
    timestamp2 = round(run_time/100);
    timestamp3 = round(run_time/10);
    timestamp4 = round(run_time);

    % snapshots in simulation time
    snap1 = (timestamp1)*tdelt;
    snap2 = (timestamp2)*tdelt;
    snap3 = (timestamp3)*tdelt;
    snap4 = (timestamp4)*tdelt;

    % times in years after applying timescale
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
        num2str(cold_temp(i)) '$^\circ$C Sink'],'FontSize',25,'Interpreter', ...
        'latex')
    xlabel('Horizontal Position','FontSize',25,'Interpreter','latex')
    ylabel('Porosity $\phi$','FontSize',25,'Interpreter','latex')
    colororder([[4/5 1 1];[3/4 4/5 0.9];[2/3 3/4 0.8];[1/2 2/3 0.7];[0 1/2 0.6]])
    legend('Initial',[num2str(time1) ' yrs after Initial'],[num2str(time2) ' yrs after Initial'],[num2str(time3) ' yrs after Initial'],[num2str(time4) ' yrs after Initial'],'fontsize',20, 'interpreter','latex','location','southeast')
    hold off

    % attempt to save files in server folders, need to ponder further
    % filename = sprintf('phi_w%0.1f_c%0.1f.mat',warm_temp(i),cold_temp(i));
    % save(['/\jumbo\ice\infiltration\supraglcial_lake_snowplug_jess\porosity_data\' filename])
    % filename = sprintf('phi_w%0.1f_c%0.1f.fig',warm_temp(i),cold_temp(i));
    % savefig(filename)
    % filename = sprintf('phi_w%0.1f_c%0.1f.png',warm_temp(i),cold_temp(i));
    % saveas(gcf,filename)
    % close all

    % save files in current folder, can save as .fig if desired
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
        num2str(cold_temp(i)) '$^\circ$C Sink'],'FontSize',25,'Interpreter', ...
        'latex')
    xlabel('Horizontal Position','FontSize',25,'Interpreter','latex')
    ylabel('Water Temperature $\theta_w$','FontSize',25,'Interpreter','latex')
    colororder([[4/5 1 1];[3/4 4/5 0.9];[2/3 3/4 0.8];[1/2 2/3 0.7];[0 1/2 0.6]])
    legend('Initial',[num2str(time1) ' mins after Initial'],[num2str(time2) ' mins after Initial'],[num2str(time3) ' mins after Initial'],[num2str(time4) ' yrs after Initial'],'fontsize',20,'interpreter','latex')
    hold off

    % save
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
        num2str(cold_temp(i)) '$^\circ$C Sink'],'FontSize',25,'Interpreter', ...
        'latex')
    xlabel('Horizontal Position','FontSize',25,'Interpreter','latex')
    ylabel('Snow Temperature $\theta_i$','FontSize',25,'Interpreter','latex')
    colororder([[4/5 1 1];[3/4 4/5 0.9];[2/3 3/4 0.8];[1/2 2/3 0.7];[0 1/2 0.6]])
    legend('Initial',[num2str(time1) ' mins after Initial'],[num2str(time2) ' mins after Initial'],[num2str(time3) ' mins after Initial'],[num2str(time4) ' yrs after Initial'],'fontsize',20,'interpreter','latex','location','southeast')
    hold off

    % save
    filename = sprintf('theta_i_w%0.1f_c%0.1f.mat',warm_temp(i),cold_temp(i));
    save(filename)
    % filename = sprintf('theta_i_w%0.1f_c%0.1f.fig',warm_temp(i),cold_temp(i));
    % savefig(filename)
    filename = sprintf('theta_i_w%0.1f_c%0.1f.png',warm_temp(i),cold_temp(i));
    saveas(gcf,filename)
    close all

    % water depth
    figure
    plot(x, h(1,:),'linewidth',2)
    hold on
    plot(x, h(timestamp1,:),'linewidth',2)
    plot(x,h(timestamp2,:),'linewidth',2)
    plot(x,h(timestamp3,:),'linewidth',2)
    plot(x, h(timestamp4,:),'linewidth',2)
    title(['Water Depth for ' num2str(warm_temp(i)) '$^\circ$C Source, ' ...
        num2str(cold_temp(i)) '$^\circ$C Sink'],'FontSize',25,'Interpreter', ...
        'latex')
    xlabel('Horizontal Position','FontSize',25,'Interpreter','latex')
    ylabel('Depth $h$','FontSize',25,'Interpreter','latex')
    colororder([[4/5 1 1];[3/4 4/5 0.9];[2/3 3/4 0.8];[1/2 2/3 0.7];[0 1/2 0.6]])
    legend('Initial',[num2str(time1) ' yrs after Initial'],[num2str(time2) ' yrs after Initial'],[num2str(time3) ' yrs after Initial'],[num2str(time4) ' yrs after Initial'],'fontsize',20,'interpreter','latex','location','southeast')
    hold off

    % save
    filename = sprintf('h_w%0.1f_c%0.1f.mat',warm_temp(i),cold_temp(i));
    save(filename)
    % filename = sprintf('h_w%0.1f_c%0.1f.fig',warm_temp(i),cold_temp(i));
    % savefig(filename)
    filename = sprintf('h_w%0.1f_c%0.1f.png',warm_temp(i),cold_temp(i));
    saveas(gcf,filename)
    close all

end