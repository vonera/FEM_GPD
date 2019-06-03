%% this script is used for generating the data for the wrr paper
load Precipitation_Externals_full_1980_2013.mat

t_start = 2;
t_end = length(precipitation_externals_full{1,3})-19*24+20;
tmp_predict(:,1) = precipitation_externals_full(:,1);

for j = 1: size(precipitation_externals_full,1)
    tmp_predict{j,2} = precipitation_externals_full{j,2}(:, t_start:t_end);
    tmp_predict{j,3} = precipitation_externals_full{j,3}(:, t_start:t_end);
    tmp_predict{j,4} = precipitation_externals_full{j,4}(:, t_start:t_end);
end


%%% normalize the local covariates to [-1 1]
% for j = 1: size(tmp_predict,1)
%     for s = 1:size(tmp_predict{j,4},1)
%         tmp_predict{j,4}(s,:) =  2*( tmp_predict{j,4}(s,:) - max( tmp_predict{j,4}(s,:) ) )./...
%             ( max( tmp_predict{j,4}(s,:) ) - min( tmp_predict{j,4}(s,:) ) ) + 1;
%     end
% end

%%% normalize the local covariates by cova ./ max(abs(cova)) 
for j = 1: size(tmp_predict,1)
    for s = 1:size(tmp_predict{j,4},1)
        tmp_predict{j,4}(s,:) =  tmp_predict{j,4}(s,:) ./ max( tmp_predict{j,4}(s,:) );
    end
end

tmp_fit(:,1) = tmp_predict(:,1);
for j = 1: size(tmp_predict,1)
    tmp_fit{j,2} = tmp_predict{j,2}(:,1:281281);
    tmp_fit{j,3} = tmp_predict{j,3}(:,1:281281);
    tmp_fit{j,4} = tmp_predict{j,4}(:,1:281281);
end



%% estimate daily accumulated precipitation
accum_precip = cell(17,1);
daily_time = cell(17,1);
cova_1 = cell(17,1);
cova_2 = cell(17,1);

for j = 1: size(tmp_predict,1)
    
    tmp = reshape(tmp_fit{j,3}(1,1:24*11720),24, 11720);
    accum_precip{j} = [sum(tmp,1) sum(tmp_fit{j,3}(24*11720+1:end))];
    daily_time{j} = tmp_fit{j,2}(1:24:end);
    
    tmp = reshape(tmp_fit{j,4}(1,1:24*11720),24, 11720);
    cova_1{j} = [sum(tmp,1)./24 sum(tmp_fit{j,4}(1,24*11720+1:end))/1];
    
    tmp = reshape(tmp_fit{j,4}(2,1:24*11720),24, 11720);
    cova_2{j} = [sum(tmp,1)./24 sum(tmp_fit{j,4}(2,24*11720+1:end))/1];
    
    tmp = reshape(tmp_fit{j,4}(3,1:24*11720),24, 11720);
    cova_3{j} = [sum(tmp,1)./24 sum(tmp_fit{j,4}(1,24*11720+1:end))/1];
    
    tmp = reshape(tmp_fit{j,4}(4,1:24*11720),24, 11720);
    cova_4{j} = [sum(tmp,1)./24 sum(tmp_fit{j,4}(2,24*11720+1:end))/1];
end


%% daily mean covariates
load('climateIndices_fit_logCo2_diffnorm')
cova = covariates_201302{2,2};
daily_mean_cova = zeros(9,11721);

for s = 1: size(cova,1)
    tmp = reshape(cova(s,1:24*11720), 24, 11720);
    daily_mean_cova(s,:) = [sum(tmp,1)./24 sum(cova(s,24*11720+1:end))/1];
end





%%
%~~~~~~~~~~~~~~~~~~~~~~~~~
QL = 0.98;
adaptSTEP = 0;
QL_eps = 0;%0.00001;
%~~~~~~~~~~~~~~~~~~~~~~~~~
%%
daily_accumulated_precipitation = cell(18,8);

daily_accumulated_precipitation{1,1} = 'city';
daily_accumulated_precipitation{1,2} = 'time_precip';
daily_accumulated_precipitation{1,3} = 'accum_precip';
daily_accumulated_precipitation{1,4} = 'daily_cova';
daily_accumulated_precipitation{1,5} = ['threshold_QL_',num2str(QL)];
daily_accumulated_precipitation{1,6} = ['time_pot_QL_',num2str(QL)];
daily_accumulated_precipitation{1,7} = ['pot_QL_',num2str(QL)];
daily_accumulated_precipitation{1,8} = ['pot_cova_QL_',num2str(QL)];
daily_accumulated_precipitation{1,9} = ['pos_switch_',num2str(adaptSTEP)];
%daily_accumulated_precipitation{1,10} = ['Btime_pot_QL_',num2str(QL)];


daily_accumulated_precipitation(2:18,1) = tmp_fit(:,1);


index = cell(17,1);
from = cell(17,1);

for j = 2:18
    daily_accumulated_precipitation{j,5} = quantile(accum_precip{j-1},QL);
    index{j-1,1} = find( accum_precip{j-1} > daily_accumulated_precipitation{j,5} + QL_eps);
    
    from{j-1,1} = 1; 
    if ~isempty( find(index{j-1,1} - 2*365 <= 0))
        from{j-1,1} = find(index{j-1,1} - 2*365 <= 0) + 1;
    end
end



for j = 2:18
    
    index_ = index{j-1,1}(from{j-1}(end):end);
    
    daily_accumulated_precipitation{j,2} = daily_time{j-1};
    daily_accumulated_precipitation{j,3} = accum_precip{j-1};
    daily_accumulated_precipitation{j,4} = [cova_1{j-1}; cova_2{j-1}; cova_3{j-1};cova_4{j-1};...
                                            daily_mean_cova;];
    
    % pot
    daily_accumulated_precipitation{j,6} = daily_time{j-1}(index_);
    daily_accumulated_precipitation{j,7} = accum_precip{j-1}(index_);
    daily_accumulated_precipitation{j,8} = [daily_accumulated_precipitation{j,4}(1:10,index_);...
                                            daily_accumulated_precipitation{j,4}(10,index_ - 90);...
                                            daily_accumulated_precipitation{j,4}(10,index_ - 365);...    
                                            daily_accumulated_precipitation{j,4}(10,index_ - 2*365);...
                                            daily_accumulated_precipitation{j,4}(11:end,index_)];
   
    
    % get adaprit position switches
    time_tmp = num2str(daily_accumulated_precipitation{j,6}(:));
    T = size(time_tmp,1);
    time_m = zeros(1,T);
    for t = 1:T
        time_m(t) = str2double(time_tmp(t,1:8));
    end
    
    diff_index = find( diff(time_m) >= adaptSTEP);
    daily_accumulated_precipitation{j,9} = unique( [1 diff_index+1 T] );
    
    %binary time
    %daily_accumulated_precipitation{j,10} = 0*daily_accumulated_precipitation{j,2};
    %daily_accumulated_precipitation{j,10}(index) = 1;
end

%%
daily_accumlated_precip_event_occurence = cell(18,2);
daily_accumlated_precip_event_occurence(2:18,1) = tmp_fit(:,1);
daily_accumlated_precip_event_occurence{1,1} = 'city';
daily_accumlated_precip_event_occurence(1,2) = daily_accumulated_precipitation(1,7);

daily_accumlated_precip_event_occurence(2:end,2) = daily_accumulated_precipitation(2:end,7);


names_covariates_fit = ['Wind','Temp','Humid','Sund', covariates_201302{1,2}(1:5),...
                        'ENSO', 'ENSO_3','ENSO_12','ENSO_24', 'CO2', 'log(CO2)', 'MJOI', 'MJOII'];

% store data
tmp = daily_accumulated_precipitation;

daily_accumulated_precipitation = []; 
daily_accumulated_precipitation.data = tmp; 
daily_accumulated_precipitation.names_cova = names_covariates_fit;
daily_accumulated_precipitation.description = ['This data set contains daily accumulated precipitation over 17',10,...
    'different locations in Switzerland. The names of locations',10,...
    'are stored in daily_accumulated_precipitation.data together with',10,...
    ['precipitation, exceedances over the ', num2str(QL),' Quantile threshold'],10,...
    'and the covariates. The names of the covariates are stored',10,...
    'in the daily_accumulated_precipitation.names_cova. All covas',10,...
    'i.e. daily_cova, do not include delayed ENSO'];

save(['daily_accumulated_precipitation_eps', num2str(QL_eps),'_QL',num2str(QL),'_delayedENSO_diffnorm.mat'], 'daily_accumulated_precipitation');


