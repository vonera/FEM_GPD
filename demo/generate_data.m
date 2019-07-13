function [Pcell, NLL_gp, Gamma, input] = generate_data()



tmp1 = ones(1,250);
tmp2 = ones(1,250);

Gamma = [tmp1 ~tmp2 tmp1 ~tmp2]; 
Gamma = [Gamma; ~Gamma]';

T = size(Gamma,1);
N = 1;

%------------------ Pcell 1 for GPD ---------------------------------------
Pcell{1,1}(1,1:3) = [ 0.3 0.0  0.0]; % param fuer ksi
Pcell{1,1}(1,4:6) = [ 2.0  2   0.5 ];  % param fuer sigma
%------------------ Pcell 2 for GPD ---------------------------------------
Pcell{2,1}(1,1:3) = [ -0.3  0.0  0.00];
Pcell{2,1}(1,4:6) = [ 1.0  0.001  -0.5];            
                            
umask_ksi = [1 2];
umask_si = [1 2];

threshold = 0*ones(N,1);

% HERE WORKING WITH FEM_Framework, load and adjust config
cfg = FEM.getdefaultCfg('MCMC','SPATIAL_BV_Gurobi');
cfg.model.subtype = 'GPD';
cfg.model.ut_combis{1} = {umask_ksi};
cfg.model.ut_combis{2} = {umask_si};
cfg.gamma.num_fe = 100;
cfg.gamma.discrete = false;

meta = Models.MCMC.GPD.Meta(T, umask_ksi, umask_si, [], [], []);

for j = 1:N
    input(j).ut = [(1/400).*(1:T); sin((2*pi)/500.*(1:T))];
    input(j).meta = meta;   
end

% generate time series wrt gpd
input = Models.MCMC.GPD.DataGen.generate_checked(Gamma, input, Pcell, threshold, T, true);     

% generate input without threshold
for j = 1:N
        input(j).xt = input(j).xt - threshold(j);
end

ksi_ut = zeros(1,T);
si_ut = zeros(1,T);

for i = 1:size(Gamma,2)
    ksi_ut = ksi_ut + Gamma(:,i)'.*( Pcell{i,1}(1,1:3)*[ones(1,T);input(j).ut] );
    si_ut  = si_ut + Gamma(:,i)'.*( Pcell{i,1}(1,4:6)*[ones(1,T);input(j).ut] ); 
end

NLL_gp = 0;
for t=1:T
    NLL_gp = NLL_gp + gplike([ksi_ut(t),si_ut(t)],input(1).xt(t));
end

%% plot
figure('Position', [10 10 500 200]) 
hist(input(1).xt, 100)
title('Histogram')

figure('Position', [10 10 500 200])
plot(input(1).xt, '*')
ylim(minmax(input(1).xt))
xlabel('t')
ylabel('xt')
title('Xt')

figure('Position', [10 10 500 200])
yyaxis left
plot(ksi_ut); 
hold on; 
yyaxis right
plot(si_ut)
xlabel('t')
ylim([0,max(max(ksi_ut, si_ut))])
legend('ksi', 'sigma')
title('Model Parameters')

end