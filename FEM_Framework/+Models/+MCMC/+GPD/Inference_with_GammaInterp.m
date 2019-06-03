%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : Inference_with_GammaInterp.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : gpd post - inference
%
% The development was guided by Prof. Dr. Illia Horenko head of the group
%    http://icsweb.inf.unisi.ch/cms/index.php/groups/20-illia-horenko.html
% USI-ICS Horenko group
%    http://icsweb.inf.unisi.ch/cms/index.php/groups/4-group-horenko.html
%--------------------------------------------------------------------------
% LICENSE
% (c) 2015, D.Igdalov, O.Kaiser USI-ICS Horenko Group, All Rights Reserved
%--------------------------------------------------------------------------
% This Software Framework may not be distributed in parts or its entirety
% without prior written agreement by O.Kaiser or D.Igdalov.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
% DEALINGS IN THE SOFTWARE.
%-------------------------------------------------------------------------
classdef Inference < handle
    %INFERENCE this class is used for inference: plott, e.g., qq, pp, roc,
    %plots for gamma but also for estimating the correlation between the
    %locations. thereby we will consider correlation, cross correlation and
    %also non-linear correlation wrt [malik2012analysis]
    
    properties
        optRes
        iHandler
        locNames
        
        Tlist
        Tlist_cs
        gammaStat
        
        timeCell
        time_all
        lon
        lat
        
        N
        K
        
        fs
        lw
    end
    
    methods
        
        %------------------------------------------------------------------
        function self = Inference(optResult, iHandler, locNames, timeCell, coord, plot_settings)
            
            self.optRes = optResult;
            self.iHandler = iHandler;
            
            % time and coordinates information
            self.timeCell = timeCell;
            self.time_all = sort(unique( cat(2,timeCell{:}) ));
            self.lon = sort(unique( cat(1,coord{:,1}) ));
            self.lat = sort(unique( cat(1,coord{:,2}) ));
            self.locNames = locNames;
            
            self.Tlist = optResult.hidden.Tlist;
            self.N = length(self.Tlist);
            self.K = optResult.hidden.K;
            
            self.Tlist_cs = cumsum([0 self.Tlist]);
            self.gammaStat = self.optRes.hidden.getPlotableStat(self.optRes.hidden.gamma);
            
            %settings for plotting the results
            if nargin > 3 && ~isempty(plot_settings)
                self.fs = plot_settings.fs;
                self.lw = plot_settings.lw;
            else
                self.fs = 15;
                self.lw = 3;
            end
        end
        
        
        %------------------------------------------------------------------
        function z_frechet = trafoFrechet(self)
            %TRAFOFRECHET tranforms the extremes to frechet distirbuted data
            
            for j = 1:self.N
                
                % Evaluate for each location: theta(t,ut)
                [ksi, sigma, ~] = self.getTheta_s(j);
                Tj = length(ksi);
                
                % Trafo
                marginal_gpd = zeros(1,Tj);
                
                for t = 1:Tj
                    if ksi(t)~=0
                        marginal_gpd(t) = 1 - (1 + ksi(t)* self.iHandler.xt{j}(t)/sigma(t))^(-1/ksi(t));
                    else
                        marginal_gpd(t) = 1 - (1 + self.iHandler.xt{j}(t)/sigma(t));
                    end
                end
                
                z_frechet{j} = -1./log(marginal_gpd);
            end
        end
        
        
        %------------------------------------------------------------------
        function plotGamma(self)
            %PLOTGAMMA plot Gamma for all locations in a subplot
            
            gamma_ = self.K*self.gammaStat;
            scrsz = get(0,'ScreenSize');
            figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
            
            for j = 1:self.N
                subplot(self.N, 1, j)
                plot(gamma_(self.Tlist_cs(j)+1:self.Tlist_cs(j+1)),'r','LineWidth', self.lw)
                ylabel(self.locNames{j}, 'FontSize', self.fs,'FontWeight','demi','Interpreter','LaTex')
                set(gca,'YTick',[1 2],'FontSize', self.fs,'FontWeight','demi')
                if j == 1
                    title('Gamma(s,t)', 'FontSize',self.fs,'FontWeight','demi','Interpreter','LaTex')
                end
            end
        end
        
        
        %------------------------------------------------------------------
        function plotTheta(self, ksiFlag, siFlag)
            % PLOTTHETA plot theta for all locations in a subplot
            
            scrsz = get(0,'ScreenSize');
            figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
            set(0, 'defaultTextInterpreter', 'LaTex');
            
            
            for j = 1:self.N
                [ksi, sigma, ~] = self.getTheta_s(j);
                
                if ksiFlag && siFlag
                    subplot(self.N, 1, j)
                    plot(ksi,'r','LineWidth', 3);
                    hold on;
                    plot(sigma,'b','LineWidth', 3)
                    ylabel(self.locNames{j}, 'FontSize', self.fs,'FontWeight','demi')
                elseif ksiFlag
                    subplot(self.N, 1, j)
                    plot(ksi,'r','LineWidth', 3);
                    ylabel(self.locNames{j}, 'FontSize',self.fs,'FontWeight','demi')
                elseif siFlag
                    subplot(self.N, 1, j)
                    plot(sigma,'b','LineWidth', 3)
                    ylabel(self.locNames{j}, 'FontSize',self.fs,'FontWeight','demi')
                end
                
                
                if j == 1
                    if ksiFlag && siFlag
                        title('Scale and Shape parameter', 'FontSize',self.fs,'FontWeight','demi')
                    elseif ksiFlag
                        title('Shape parameter', 'FontSize',self.fs,'FontWeight','demi')
                    elseif siFlag
                        title('Scale parameter: $\sigma(t)$', 'FontSize',self.fs,'FontWeight','demi')
                    end
                end
            end
        end
        
        
        %------------------------------------------------------------------
        function [t_T1, yt] = plotQP(self, ppFlag, qqFlag, threshold_all)
            %PLOTQP provides the QQ and the PP plot for the fitted gpd
            
            t_T1 = cell(1,self.N);
            yt = cell(1,self.N);
            
            % Trafo to standard exp data
            for j = 1:self.N
                [ksi, sigma, ~] = self.getTheta_s(j);
                yt{j} = (1./ksi) .* log( 1 + ksi .* (self.iHandler.xt{j} - threshold_all(j))./sigma );
                yt{j} = sort(yt{j});
                t = 1:length(ksi);
                t_T1{j} = t./(length(ksi) +1);
            end
            
            % PP plot
            if ppFlag
                scrsz = get(0,'ScreenSize');
                figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
                set(0, 'defaultTextInterpreter', 'LaTex');
                for j = 1:self.N
                    subplot(ceil(self.N/5), 5, j)
                    plot( t_T1{j},  t_T1{j}, 'k'); hold on;
                    plot( t_T1{j}, 1 - exp( -yt{j} ),'r','LineWidth', 3);
                    xlabel('Empirical', 'FontSize',self.fs,'FontWeight','demi')
                    ylabel('Observed', 'FontSize',self.fs,'FontWeight','demi')
                    title(['PP plot for ',self.locNames{j}], 'FontSize',self.fs,'FontWeight','demi')
                end
            end
            
            %QQ plot
            if qqFlag
                scrsz = get(0,'ScreenSize');
                figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
                set(0, 'defaultTextInterpreter', 'LaTex');
                for j = 1:self.N
                    subplot(ceil(self.N/5), 5, j)
                    plot( -log( 1 - t_T1{j}),-log( 1 - t_T1{j}), 'k'); hold on;
                    plot( -log( 1 - t_T1{j}), yt{j} ,'r','LineWidth', 3);
                    xlabel('Empirical', 'FontSize',self.fs,'FontWeight','demi')
                    ylabel('Observed', 'FontSize',self.fs,'FontWeight','demi')
                    title(['QQ plot for ',self.locNames{j}], 'FontSize',self.fs,'FontWeight','demi')
                end
            end
        end
        
        
        %------------------------------------------------------------------
        function plotRocCdf(self, event, threshold_all)
            %PLOTROCCDF evaluates the probability for an event and compares
            %it to the empirical probability.
            
            emp_cdf = cell(1,self.N);
            cdf = cell(1,self.N);
            
            for j = 1:self.N
                Tj = self.Tlist(j);
                tmp_cdf = zeros(1,Tj);
                emp_cdf{j} = zeros(1,Tj);
                
                for i = 1:self.K
                    Gammaj = self.optRes.hidden.gamma(self.Tlist_cs(j)+1:self.Tlist_cs(j+1),i)';
                    ksi = self.optRes.param.theta(i,self.optRes.param.ksi_mask)* [ones(1,Tj); self.iHandler.ut{j}(self.optRes.param.meta.umask_ksi,:)];
                    sigma = self.optRes.param.theta(i,self.optRes.param.si_mask)* [ones(1,Tj); self.iHandler.ut{j}(self.optRes.param.meta.umask_si,:)];
                    
                    % compute probabilities
                    value(ksi~=0) =  ksi(ksi~=0).*( event(j)- threshold_all(j) ) ./ sigma(ksi~=0);
                    tmp_cdf(ksi~=0) = tmp_cdf(ksi~=0) + Gammaj(ksi~=0).*(1 + value(ksi~=0)) .^ (-1./ksi(ksi~=0));
                end
                
                cdf{j} =  tmp_cdf;
                emp_cdf{j}(self.iHandler.xt{j} > event(j)) = 1;
            end
            
            scrsz = get(0,'ScreenSize');
            figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
            set(0, 'defaultTextInterpreter', 'LaTex');
            
            for j = 1:self.N
                [tpr,fpr,~] =  roc(emp_cdf{j},cdf{j});
             
                subplot(ceil(self.N/5), 5, j)
                plot(fpr,tpr,'b','LineWidth',self.lw);
                xlabel('FPR', 'FontSize',self.fs,'FontWeight','demi')
                ylabel('TPR', 'FontSize',self.fs,'FontWeight','demi')
                title(['P[$X >$ ',num2str(event(j)),']',' for ',self.locNames{j}], 'FontSize',self.fs,'FontWeight','demi')
            end
        end
        
        
        %------------------------------------------------------------------
        function returnLevel = computeReturnLevel(self, threshold_all, p, length_origXt)
            %COMPUTERETURNLEVEL provides the effective returnlevel for the
            %return period p (effective means nonstationary)
            
            returnLevel = cell(1,self.N);
            
            for j=1:self.N
                psi_u = self.Tlist(j) / length_origXt;
                Tj = self.Tlist(j);
                [ksi, sigma, ~] = self.getTheta_s(j);
                
                for l=1:length(p)
                    for t = 1:Tj
                        if ksi(t)~=0
                            returnLevel{j}(l,t) = threshold_all(j) + sigma(t)/ksi(t)*( ( p(l)*psi_u )^ksi(t) - 1 );
                        else
                            returnLevel{j}(l,t) = threshold_all(j) + sigma(t)*log( p(l)*psi_u );
                        end
                    end
                end
            end
        end
        
        
        %------------------------------------------------------------------
        function [ksi, sigma, Gammaj] = getTheta_s(self, j)
            %GETTHETA_S get for location s parameter ksi and sigma and Gamma
            
            Tj = self.Tlist(j);
            Gammaj = zeros(self.K,Tj);
            sigma = zeros(1,Tj);
            ksi = zeros(1,Tj);
            
            for i = 1:self.K
                Gammaj(i,:) = self.optRes.hidden.gamma(self.Tlist_cs(j)+1:self.Tlist_cs(j+1),i)';
                ksi = ksi + Gammaj(i,:).*(self.optRes.param.theta(i,self.optRes.param.ksi_mask)* [ones(1,Tj); self.iHandler.ut{j}(self.optRes.param.meta.umask_ksi,:)]);
                sigma = sigma + Gammaj(i,:).*(self.optRes.param.theta(i,self.optRes.param.si_mask)* [ones(1,Tj); self.iHandler.ut{j}(self.optRes.param.meta.umask_si,:)]);
            end
        end
        
        
        %------------------------------------------------------------------
        function [GammaFull, lon_, lat_] = interpolateGamma(self, coord, maxIter, gridS)
            %INTEROPLATEGAMMA is used to interpolate the Gamma between the
            %the coordinates of the corresponding locations. The interpolated
            %region is a squered grid. The interpolation is not max-stable.
            %The interpolation is based on smoothing splines.
            
            if ~isempty(gridS)
                [lon_, lat_] = self.fineGrid(gridS, self.lon, self.lat);
            else
                lon_ = self.lon;
                lat_ = self.lat;
            end
            
            Gamma = NaN * ones(length(lon_), length(lat_), length(self.time_all));
            gamma_ = self.K * self.gammaStat;
            
            for j = 1:self.N
                time_j = self.timeCell{j};
                index_lon_j = find(coord{j,1} == lon_);
                index_lat_j = find(coord{j,2} == lat_);
                Gammaj = gamma_(self.Tlist_cs(j)+1:self.Tlist_cs(j+1));
                
                for t = 1:length(Gammaj)
                    index_jt = find(time_j(t) == self.time_all);
                    Gamma(index_lon_j, index_lat_j, index_jt) = Gammaj(t); %#ok<FNDSB>
                end
            end
            
            % use interpolation inpaintn.m
            if nargin > 3 && ~isempty(maxIter)
                numIter = maxIter;
            else
                numIter = 100;
            end
            
            GammaFull = Models.MCMC.GPD.Tools.inpaintn(Gamma,numIter);
        end
        
        
        %------------------------------------------------------------------
        function corrcoefGamma = computeCorrCoef(self, equidFlag, surfFlag)
            
            if equidFlag
                aGamma = self.allignGammaEqui();
            else
                aGamma = self.allignGamma();
            end
            corrcoefGamma = corrcoef(aGamma');
            
            % plot
            if surfFlag
                scrsz = get(0,'ScreenSize');
                figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
                set(0, 'defaultTextInterpreter', 'LaTex');
                surf(corrcoefGamma);
                title('Corrcoef for Gamma','FontSize',self.fs,'FontWeight','demi')
                xlim([1 self.N])
                ylim([1 self.N])
                set(gca, 'XTick', 1:self.N, 'XTickLabel', self.locNames(1:self.N), 'FontSize',self.fs,'FontWeight','demi')
                set(gca, 'YTick', 1:self.N, 'YTickLabel', self.locNames(1:self.N),'FontSize',self.fs,'FontWeight','demi')
                view([-90.5 90]);
                colorbar('location', 'WestOutside','FontSize',self.fs,'FontWeight','demi');
                caxis([0 1])
            end
        end
        
        
        %------------------------------------------------------------------
        function corrcoefGamma = computeCorrCoefCl(self, cl, equidFlag, surfFlag)
            
            if equidFlag
                aGamma = self.allignGammaEqui();
            else
                aGamma = self.allignGamma();
            end
            
            index = 1:self.K;
            index_cl = find(index ~= cl);
            
            for k = index_cl
                aGamma( aGamma == k ) = 0;
            end
            
            corrcoefGamma = corrcoef(aGamma');
            
            % plot if desired
            if surfFlag
                scrsz = get(0,'ScreenSize');
                figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
                set(0, 'defaultTextInterpreter', 'LaTex');
                surf(corrcoefGamma);
                title(['Corrcoef for Gamma, cluster',num2str(cl)],'FontSize',self.fs,'FontWeight','demi')
                xlim([1 self.N])
                ylim([1 self.N])
                set(gca, 'XTick', 1:self.N, 'XTickLabel', self.locNames(1:self.N), 'FontSize',self.fs,'FontWeight','demi')
                set(gca, 'YTick', 1:self.N, 'YTickLabel', self.locNames(1:self.N),'FontSize',self.fs,'FontWeight','demi')
                view([-90.5 90]);
                colorbar('location', 'WestOutside','FontSize',self.fs,'FontWeight','demi');
                caxis([0 1])
            end
        end
        
        
        %------------------------------------------------------------------
        function [xcorrGamma, lags] = computeXCorr(self, equidFlag, timeLag, stemFlag, sliceFlag)
            
            if equidFlag
                aGamma = self.allignGammaEqui();
            else
                aGamma = self.allignGamma();
            end
            [xcorrGamma, lags] = xcorr(aGamma', timeLag, 'coef');
            
            % plot
            if stemFlag
                figure;
                stem(lags, xcorrGamma);
            end
            
            if sliceFlag
                xcorr_slice = zeros(self.N, self.N, size(xcorrGamma,1));
                
                for t = 1:size(xcorrGamma,1)
                    xcorr_slice(:,:,t) = reshape( squeeze(xcorrGamma(t,:)), 17, 17 );
                end
                
                scrsz = get(0,'ScreenSize');
                figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
                set(0, 'defaultTextInterpreter', 'LaTex');
                h = slice(xcorr_slice, [], [], 1:size(xcorr_slice,3));
                set(h,'EdgeColor','none','FaceColor','interp')
                alpha(.3)
            end
        end
    end
    
    
    methods (Access = 'protected')
        
        
        %------------------------------------------------------------------
        function [coord_lon_f, coord_lat_f] = fineGrid(~, gridS, lon, lat)
            
            % refine the lon
            coord_lon_f = lon(1);
            diff_crd_lon = diff(lon);
            
            for t = 1:length(diff_crd_lon)
                if diff_crd_lon(t) > gridS
                    n = round(round((lon(t+1)-lon(t))/gridS));
                    step_size = ( lon(t+1)-lon(t) )/n;
                    coord_lon_f = [coord_lon_f lon(t):step_size:lon(t+1)];
                end
                coord_lon_f = [coord_lon_f lon(t)];
            end
            coord_lon_f = unique([coord_lon_f lon(t+1)]);
            
            % refine the lat
            coord_lat_f = lat(1);
            diff_crd_lat = diff(lat);
            
            for t = 1:length(diff_crd_lat)
                if diff_crd_lat(t) > gridS
                    n = round(round((lat(t+1)-lat(t))/gridS));
                    step_size = ( lat(t+1)-lat(t) )/n;
                    coord_lat_f = [coord_lat_f lat(t):step_size:lat(t+1)]; %#ok<*AGROW>
                end
                coord_lat_f = [coord_lat_f lat(t)];
            end
            coord_lat_f = unique([coord_lat_f lat(t+1)]);
        end
        
        
        %------------------------------------------------------------------
        function aGamma = allignGamma(self)
            
            aGamma = zeros(self.N, length(self.time_all));
            
            for j = 1:self.N
                time_j = self.timeCell{j};
                gamma_j = self.K*self.gammaStat(self.Tlist_cs(j)+1:self.Tlist_cs(j+1));
                
                for t = 1:length(gamma_j)
                    aGamma(j, time_j(t) == self.time_all) = gamma_j(t);
                end
            end
        end
        
        
        %------------------------------------------------------------------
        function aGammaEquid = allignGammaEqui(self)
            
            timeEqui = self.time_all(1):self.time_all(end);
            aGammaEquid = zeros(self.N, length(timeEqui));
            
            for j = 1:self.N
                time_j = self.timeCell{j};
                gamma_j = self.K*self.gammaStat(self.Tlist_cs(j)+1:self.Tlist_cs(j+1));
                
                for t = 1:length(gamma_j)
                    aGammaEquid(j, time_j(t) == self.time_all) = gamma_j(t);
                end
            end
        end
    end
    
    
end

