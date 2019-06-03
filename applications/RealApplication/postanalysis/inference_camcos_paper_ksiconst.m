%for QL=0.98
clear all
close all

S_crd = load('../coordinates');
S =load('../daily_accumulated_precipitation_0.98_delayedENSO_diffnorm.mat');
load('../ksi_const_potSwiss_QL0.98_lasso0_AICc24587.1269.mat')

%% use infHandler for post analysis
coord = S_crd.coordinates(2:end,4:5);
timeCell = S.daily_accumulated_precipitation.data(2:end,6);
locNames = S_crd.coordinates(2:end,1);

infHandler = Models.MCMC.GPD.Inference(optRes, iHandler, locNames, timeCell, coord, []);


%% plot theta for location Lugano, Bern, Alt
plotParams(infHandler, timeCell, {'BER','BAS','LUG','ALT'})


%% get confidence intervalls for theta
conf_int = infHandler.getConf4Theta();

%% correlation & synchronisation clusterwise

[~, strengthES_cl1] = infHandler.computeES(1);
[~, strengthES_cl2] = infHandler.computeES(2);

% plot the synchronisation measure
myGCA.fsize = 24;

% plot sorted synchronisation measure
A_coord = cell2mat(S_crd.coordinates(2:end,2:3));
[B,index] = sortrows(A_coord,[2,1]);

index = fliplr(index')';

plotM(strengthES_cl1(index,index), locNames(index), myGCA)
plotM(strengthES_cl2(index,index), locNames(index), myGCA)


%% plot gamma for locations
% 'GUT',    'BAS',    'KLO',    'FAH',    'GLA'
% 'NAP',    'NEU',    'BER',    'ALT',    'CHU',    'DAV'
% 'SCU',    'ULR',    'CGI',    'GVE',    'SIO',    'LUG'

gamma4_loc = {'GUT',    'BAS',    'KLO',    'FAH',    'GLA'};
Ntmp = length(gamma4_loc);

figure;
for j = 1:Ntmp    
    index_loc =  find(ismember(infHandler.locNames,gamma4_loc{j}));
    [~, ~, Gammaj] = infHandler.getTheta_s(index_loc);
    Gammaj = Gammaj + 1; 
    subplot(Ntmp, 1, j)
    plot(infHandler.serialDate_nr{index_loc}, Gammaj(2,:),':.k','MarkerSize', 8*infHandler.lw)
    datetick('x','yy','keeplimits');
    ylim([0.8 2.2])
    ylabel(infHandler.locNames{index_loc}, 'FontSize', infHandler.fs,'FontWeight','demi','Interpreter','LaTex')
    set(gca,'YTick',[1 2],'FontSize', infHandler.fs,'FontWeight','demi')
    if j == 1
        title('Gamma(s,t)', 'FontSize',infHandler.fs,'FontWeight','demi','Interpreter','LaTex')
    end
end





