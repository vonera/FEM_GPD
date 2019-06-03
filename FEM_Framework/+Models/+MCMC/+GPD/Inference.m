%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : Inference_with_GammaInterp.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : gpd post inference
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
        serialDate_nr
        lon
        lat
        
        N
        NN
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
            self.NN = self.N^2;
            self.K = optResult.hidden.K;
            
            self.Tlist_cs = cumsum([0 self.Tlist]);
            self.gammaStat = self.optRes.hidden.getPlotableStat(self.optRes.hidden.gamma);
            self.serialDate_nr = self.getTimeEx();
            
            %settings for plotting the results
            if nargin > 3 && ~isempty(plot_settings)
                self.fs = plot_settings.fs;
                self.lw = plot_settings.lw;
            else
                self.fs = 20;
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
                plot(self.serialDate_nr{j}, gamma_(self.Tlist_cs(j)+1:self.Tlist_cs(j+1)),'.k','MarkerSize', 3*self.lw)
                datetick('x','yy','keeplimits');
                ylim([0.8 2.2])
                ylabel(self.locNames{j}, 'FontSize', self.fs,'FontWeight','demi','Interpreter','LaTex')
                set(gca,'YTick',[1 2],'FontSize', self.fs,'FontWeight','demi')
                if j == 1
                    title('Gamma(s,t)', 'FontSize',self.fs,'FontWeight','demi','Interpreter','LaTex')
                end
            end
        end
        
        
        %------------------------------------------------------------------
        function plotTheta(self, ksiFlag, siFlag, ksiOrig, siOrig)
            % PLOTTHETA plot theta for all locations in a subplot
            
            scrsz = get(0,'ScreenSize');
            figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
            set(0, 'defaultTextInterpreter', 'LaTex');
            
            
            for j = 1:self.N
                [ksi, sigma, ~] = self.getTheta_s(j);
                
                if ksiFlag && siFlag
                    subplot(self.N, 1, j)
                    % plot ksi
                    plot(self.serialDate_nr{j},ksi,'r','LineWidth', 3);
                    datetick('x','yy','keeplimits');
                    hold on;
                    % plot sigma
                    plot(self.serialDate_nr{j},sigma,'b','LineWidth', 3)
                    datetick('x','yy','keeplimits');
                    ylabel(self.locNames{j}, 'FontSize', self.fs,'FontWeight','demi')
                    set(gca,'FontSize', self.fs,'FontWeight','demi')
                elseif ksiFlag
                    subplot(self.N, 1, j)
                    plot(self.serialDate_nr{j},ksi,'r','LineWidth', 3);
                    datetick('x','yy','keeplimits');
                    hold on;
                    if nargin >= 4
                       plot(self.serialDate_nr{j}, ksiOrig{j},'--b','LineWidth', 2); 
                       datetick('x','yy','keeplimits');
                    end
                    ylabel(self.locNames{j}, 'FontSize',self.fs,'FontWeight','demi')
                    set(gca,'FontSize', self.fs,'FontWeight','demi')
                elseif siFlag
                    subplot(self.N, 1, j)
                    if nargin == 5
                        plot(siOrig{j},'-.r','LineWidth', 3);%plot(self.serialDate_nr{j}, siOrig{j},'--r','LineWidth', 3);
                        %datetick('x','yy','keeplimits');
                    end
                    hold on;
                    plot(sigma,'-k','LineWidth', 3);%plot(self.serialDate_nr{j}, sigma,'-k','LineWidth', 3)
                    %datetick('x','yy','keeplimits');
                    ylabel(self.locNames{j}, 'FontSize',self.fs,'FontWeight','demi')
                    set(gca,'FontSize', self.fs,'FontWeight','demi')
                    xlim([0,565])
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
            xlabel('time', 'FontSize',self.fs,'FontWeight','demi')
        end
        
        
        %------------------------------------------------------------------
        function [t_T1, yt] = plotQP(self, ppFlag, qqFlag, threshold_all, QQ_bands)
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
                    plot(t_T1{j},  t_T1{j}, 'k','LineWidth', 3); 
                    datetick('x','yy','keeplimits');
                    hold on;
                    plot(t_T1{j}, 1 - exp( -yt{j} ),'*r','LineWidth', 3);
                    datetick('x','yy','keeplimits');
                    xlabel('Empirical', 'FontSize',self.fs,'FontWeight','demi')
                    ylabel('Model', 'FontSize',self.fs,'FontWeight','demi')
                    title(['PP plot for ',self.locNames{j}], 'FontSize',self.fs,'FontWeight','demi')
                end
            end
            
            %QQ plot
            if qqFlag
                if isempty(QQ_bands)
                    scrsz = get(0,'ScreenSize');
                    figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
                    set(0, 'defaultTextInterpreter', 'LaTex');
                    for j = 1:self.N
                        subplot(ceil(self.N/5), 5, j)
                        plot( yt{j},yt{j}, 'k','LineWidth', 3); hold on;
                        plot( yt{j}, -log( 1 - t_T1{j}),'*r','LineWidth', 3);
                        xlabel('Model', 'FontSize',self.fs,'FontWeight','demi')
                        ylabel('Empirical', 'FontSize',self.fs,'FontWeight','demi')
                        title(['QQ plot for ',self.locNames{j}], 'FontSize',self.fs,'FontWeight','demi')
                    end
                else
                    grey=[0.4,0.4,0.4];
                    scrsz = get(0,'ScreenSize');
                    figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
                    set(0, 'defaultTextInterpreter', 'LaTex');
                    for j = 1:self.N
                        subplot(ceil(self.N/5), 5, j)
                        plot( yt{j},yt{j}, 'k','LineWidth', 3); hold on;
                        plot( yt{j}, -log( 1 - t_T1{j}),'.r','MarkerSize',22);
                        plot(QQ_bands{j}(:,1),-log( 1 - t_T1{j}), '--','Color',grey, 'LineWidth', 2);
                        plot(QQ_bands{j}(:,2),-log( 1 - t_T1{j}), '--','Color',grey, 'LineWidth', 2);
                        xlabel('Model', 'FontSize',self.fs,'FontWeight','demi')
                        ylabel('Empirical', 'FontSize',self.fs,'FontWeight','demi')
                        title(['QQ plot for ',self.locNames{j}], 'FontSize',self.fs,'FontWeight','demi')
                    end
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
        function corrcoefGamma = computeCorrCoef(self, surfFlag, equidFlag)
            %COMPUTECORRCOEF estimates the relationship of the affiliation
            % vetor among locations. For this we have two possibilities.
            % First: the affiliation vectors are first aligned to real time line,
            % such that A(s,t) = 0,1,2,..,K. Where 0 is assigned to time steps
            % where no extreme occured. Second: affiliation vector is assigned
            % to the ocurrence of extreme events. The first approach seems to
            % be more intuitive, and thus is the default.
            
            if equidFlag
                aGamma = self.allignGammaEqui();
            else
                aGamma = self.allignGamma();
            end
            corrcoefGamma = corrcoef(aGamma');
            
            %plot
            if surfFlag
                titleStr = 'Correlation of Affiliation vector A(s,t), wrt to corrcoef()';
                self.plotCorrMatrix(corrcoefGamma,titleStr)
            end
        end
        
        
        %------------------------------------------------------------------
        function corrcoefGamma = computeCorrCoefCl(self, cl, equidFlag, surfFlag)
            %COMPUTECORRCOEFCL 
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
            
            % plot
            if surfFlag
                titleStr = 'Correlation of Affiliation vector A(s,t), wrt to corrcoef(), clusterwise';
                self.plotCorrMatrix(corrcoefGamma,titleStr)
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
        
        
        %------------------------------------------------------------------
        function [delay_ES, strength_ES] = computeES(self, cl)
            %COMPUTEES estimates the event synchronization matrices
            %The results are taken from:
            %https://www.pik-potsdam.de/members/kurths/publikationen/2012/malik%20clim%20dyn39.pdf
            %GETDELAYES estiamtes the event synch strength/delay, matrices
            
            % input: cl for wich cluster the ES should be obtained
            % if cl = 0, then it computes the stationary ES
            if cl ~= 0
                tES = self.getClOccurence(self.serialDate_nr, cl);
            else
                tES = self.serialDate_nr;
            end
            matrixES = self.getES(tES);
            
            strength_ES = zeros(self.N,self.N);
            delay_ES = zeros(self.N,self.N);
            
            for i = 1:self.N
                for j = 1:self.N
                    nomr_const = sqrt( length(tES{i}) * length(tES{j}) );
                    strength_ES(i,j) = ( matrixES(i,j) + matrixES(j,i) ) / nomr_const;
                    delay_ES(i,j) = ( matrixES(i,j) - matrixES(j,i) ) / nomr_const;
                end
            end
        end
        
        
        %------------------------------------------------------------------
        function std_error = getConf4Theta(s)
            % estimate the confidence intervals by evaluating the numerical
            % Fischer information matrix using the matlab fucntion mlecov
            
            Ust_all = cat(2, s.iHandler.ut{:});
            Xst_all = cat(2, s.iHandler.xt{:});
            cfg = s.iHandler.cfg;
            
            std_error = zeros(s.optRes.hidden.K,length(s.optRes.param.theta));
            
            for cl=1:s.optRes.hidden.K
                
                Gamma_cl = s.optRes.hidden.gamma(:,cl);
                theta_cl = s.optRes.param.theta(cl,:);
                
                Ust_cl = Ust_all(cfg.model.ut_combis{1}{1},Gamma_cl==1);
                Xst_cl = Xst_all(Gamma_cl==1);
                
                f = @(params,data, Ust, ksi_mask,si_mask,cens,freq)s.getLL(data, params, Ust_cl, s.optRes.param.ksi_mask, s.optRes.param.si_mask);
                
                acov_cl = mlecov(theta_cl, Xst_cl, 'nloglf', f);
                
                std_error(cl,:) = sqrt(diag(acov_cl));
                
            end
        end
        
        
    end
    
    
    methods (Access = 'protected')
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
        
        
        %------------------------------------------------------------------
        function plotCorrMatrix(self, corrMat, titleStr)
            scrsz = get(0,'ScreenSize');
            figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)],...
                   'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.100000001490116 1;0 0.200000002980232 1;0 0.300000011920929 1;0 0.400000005960464 1;0 0.5 1;0 0.600000023841858 1;0 0.699999988079071 1;0 0.800000011920929 1;0 0.899999976158142 1;0 1 1;0.0909090936183929 1 0.909090936183929;0.181818187236786 1 0.818181812763214;0.272727280855179 1 0.727272748947144;0.363636374473572 1 0.636363625526428;0.454545468091965 1 0.545454561710358;0.545454561710358 1 0.454545468091965;0.636363625526428 1 0.363636374473572;0.727272748947144 1 0.272727280855179;0.818181812763214 1 0.181818187236786;0.909090936183929 1 0.0909090936183929;1 1 0;1 0.899999976158142 0;1 0.800000011920929 0;1 0.699999988079071 0;1 0.600000023841858 0;1 0.5 0;1 0.400000005960464 0;1 0.300000011920929 0;1 0.200000002980232 0;1 0.100000001490116 0;1 0 0;0.980000019073486 0 0;0.959999978542328 0 0;0.939999997615814 0 0;0.920000016689301 0 0;0.899999976158142 0 0;0.879999995231628 0 0;0.860000014305115 0 0;0.839999973773956 0 0;0.819999992847443 0 0;0.800000011920929 0 0;0.779999971389771 0 0;0.759999990463257 0 0;0.740000009536743 0 0;0.720000028610229 0 0;0.699999988079071 0 0;0.680000007152557 0 0;0.660000026226044 0 0;0.639999985694885 0 0;0.620000004768372 0 0;0.600000023841858 0 0;0.579999983310699 0 0;0.560000002384186 0 0;0.540000021457672 0 0;0.519999980926514 0 0;0.5 0 0]);
            set(0, 'defaultTextInterpreter', 'LaTex');
            imagesc(corrMat);                        
            textStrings = num2str(corrMat(:),'%0.2f');    %# Create strings from the matrix values
            textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
            [x,y] = meshgrid(1:self.N);                   %# Create x and y coordinates for the strings
            hStrings = text(x(:),y(:),textStrings(:),...  %# Plot the strings
                'HorizontalAlignment','center','FontSize',...
                 self.fs,'FontWeight','demi');
            midValue = mean(get(gca,'CLim'));             %# Get the middle value of the color range
            textCol = repmat(corrMat(:) > midValue,1,3);  %# Choose white or black for the
            
            %   text color of the strings so they can be easily seen over the background color
            set(hStrings,{'Color'},num2cell(textCol,2));  %# Change the text colors
            set(gca,'XTick',1:self.N,...                  %# Change the axes tick marks
                'XTickLabel',self.locNames,...            %#   and tick labels
                'YTick',1:self.N,...
                'YTickLabel',self.locNames,...
                'TickLength',[0 0],'FontSize',self.fs,'FontWeight','demi');
            title(titleStr, 'FontSize',self.fs,'FontWeight','demi')
        end
        
        
        %------------------------------------------------------------------
        function tES = getClOccurence(self,time_x, cl)
            % input: time_x - 1xN cell with occurence (time) of extremes
            tES = cell(1, self.N);
             
            for j = 1:self.N;
                [~, ~, gammaj] = getTheta_s(self, j);
                plotGammaj = self.K * self.optRes.hidden.getPlotableStat(gammaj');
                tES{j} = time_x{j}(plotGammaj == cl);
            end
        end
        
        
        %------------------------------------------------------------------
        function event_synch = getES(s, tES)
            %GETES gets the pairwise synronisation matrix
            event_synch = zeros(s.N, s.N);
            for i = 1:s.N
                lj = length(tES{i});
                diff_i = [inf diff(tES{i})' inf];
                
                for j = 1:s.N
                    tmp = 0;
                    
                    li = length(tES{j});
                    diff_j = [inf diff(tES{j})' inf];
                    
                    for ti = 1:lj-1
                        
                        for tj= 1:li-1
                            tlag_ij = min([diff_i(ti:ti+1), diff_j(tj:tj+1)]);
                            
                            if tES{i}(ti) == tES{j}(tj)
                                tmp = tmp + 0.5;
                            elseif ( 0 < tES{i}(ti) - tES{j}(tj) ) && ...
                                    ( tES{i}(ti) - tES{j}(tj) < tlag_ij/2 )
                                tmp = tmp + 1;
                            end
                        end
                    end
                    
                    event_synch(i,j) = tmp;
                end
            end
        end
        
        
        %------------------------------------------------------------------
        function time_ex = getTimeEx(self)
        % this function transofrms data format into matlab format.
            time_ex = cell(1,self.N);
            
            for j = 1:self.N;
                
                tmp = zeros(length(self.timeCell{j}),1);
                
                for t = 1:length(self.timeCell{j})
                    a = num2str(self.timeCell{j}(t));
                    tmp(t) = datenum([a(1:4),'-',a(5:6),'-',a(7:8)]);
                end
                
                time_ex{j} = tmp;
            end
        end
        
        
        %------------------------------------------------------------------
        function nlogpdf = getLL(~, xt_k, params, Ust_cl, ksi_mask, si_mask)

            ksi_k = params(ksi_mask);
            si_k = params(si_mask);
            
            nlogpdf = 0;
            
            for t = 1:length(xt_k)
                
                ksi = ksi_k * [1; Ust_cl(:,t)];
                si  = si_k * [1; Ust_cl(:,t)];
                
                if ksi ~= 0
                    zt = 1 + ksi * xt_k(t) / si;
                    if zt > 0
                        nlogpdf = nlogpdf + log(si) + (1+1/ksi)*log(zt);
                    end
                else
                    nlogpdf = log(si) + xt_k(t) / si;
                end
            end
        end
        
    end
    
    
end

