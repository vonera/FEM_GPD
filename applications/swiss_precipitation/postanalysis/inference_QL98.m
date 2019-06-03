%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This scipr is used for post-analysis of the results and for ploting
% figures in the paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

%% load results and data
S_crd = load('../data/coordinates');
S =load('../data/daily_accumulated_precipitation_0.98_delayedENSO_diffnorm.mat');
load('../results/potSwiss_QL0.98_lasso0_AICc24625.3443.mat')

%% initialize inference Handler (infHandler) for post-analysis
coord = S_crd.coordinates(2:end,4:5);
timeCell = S.daily_accumulated_precipitation.data(2:end,6);
locNames = S_crd.coordinates(2:end,1);

infHandler = Models.MCMC.GPD.Inference(optRes, iHandler, locNames, timeCell, coord, []);
z_frechet = infHandler.trafoFrechet();

%% get confidence intervalls for theta -> Table 1
conf_int = infHandler.getConf4Theta();


%% plot gamma for locations -> Figure 1 
gamma4_loc = {'BAS', 'DAV', 'CGI','LUG'};
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


%% plot theta for locations -> Figure 2
plotParams(infHandler, timeCell, {'BAS', 'DAV', 'CGI', 'LUG'})


%% load and plot confidence intervals for QQ plots -> Figure 3 
path = pwd;
cd frechetData_4Rstudio_QL98/qq_bands

% get the bounds for qq in a cell 17x1
cell_all_bands = cell(17,1);
for j=1:17
    filename = [num2str(j),'ycl_all','_bands.mat'];
    Sbands = load(filename);
    cell_all_bands{j} = [Sbands.qq_bands_lo Sbands.qq_bands_up];
end
cd(path)

threshold_all = 0*cat(1, S.daily_accumulated_precipitation.data{2:end,5});
infHandler.plotQP(false, true, threshold_all, cell_all_bands);


%% Figure 4 & 5 & 6 were provided by Olivia R.


%% synchronisation stationary and clusterwise -> Figure 7
corrcoefGamma_cl1 = infHandler.computeCorrCoefCl(1,true, true);
corrcoefGamma_cl2 = infHandler.computeCorrCoefCl(2,true, true);

[delayES_cl1, strengthES_cl1] = infHandler.computeES(1);
[delayES_cl2, strengthES_cl2] = infHandler.computeES(2);
[delayES, strengthES] = infHandler.computeES(0);


% plot sorted synchronisation measure
A_coord = cell2mat(S_crd.coordinates(2:end,2:3));
[B,index] = sortrows(A_coord,[2,1]);

index = fliplr(index')';

myGCA.fsize = 24;
plotM(strengthES_cl1(index,index), locNames(index), myGCA)
plotM(strengthES_cl2(index,index), locNames(index), myGCA)
plotM(strengthES(index,index), locNames(index), myGCA)

