function plotParams(infHandler, timeCell, locNames)
% this function plots the parameters for locations which are listed in the
% input parameter locNames


fsize = 24;
ms = 30;
lw = 20;
grey = [0.4,0.4,0.4];
figure;
warning('here plots for K=2')

Nloc = length(locNames);
index_j = 1;

for j = 1:Nloc

    % find index of the location
    index_loc =  find(ismember(infHandler.locNames,locNames{j}));
    index_loc = index_loc;

    % make subplot for the location
    index_time = zeros(1,length(timeCell{index_loc}));
    for t = 1:length(timeCell{index_loc})
        a = num2str(timeCell{index_loc}(t));
        index_time(t) = datenum([a(1:4),'-',a(5:6),'-',a(7:8)]);
    end
    
    [ksi_fem,si_fem,Gammaj] =infHandler.getTheta_s(index_loc);
    
    subplot(Nloc,2,index_j)
    plot(index_time(Gammaj(1,:)==1),si_fem(Gammaj(1,:)==1),'.k','MarkerSize',ms,'LineWidth',lw);hold on;
    plot(index_time(Gammaj(2,:)==1),si_fem(Gammaj(2,:)==1),'.r','MarkerSize',ms,'LineWidth',lw);
    xlim(minmax(index_time))
    datetick('x','yy','keeplimits');
    xlabel('Year','FontSize',fsize,'FontWeight','demi','Interpreter','LaTex')
    title(['Scale Parameter $\sigma(U_t)$ for location ',infHandler.locNames{index_loc}],'Interpreter','LaTex','FontSize',fsize,'FontWeight','demi')
    h = legend('$\sigma_1(U_t)$','$\sigma_2(U_t)$', 'Orientation','horizontal');
    set(h,'Interpreter','LaTex','Box','on','FontSize',fsize,'FontWeight','demi')
    ylim([0 65])
    set(gca,'YTick',[ 10  30 ],'YTickLabel',[10  30 ],'FontSize',fsize);%,'FontName','cmr12')
    index_j = index_j + 1; 
     
     
    subplot(Nloc,2,index_j)
    plot(minmax(index_time),[0 0],'-','LineWidth',1,'Color',grey);hold on;
    plot(index_time(Gammaj(1,:)==1),ksi_fem(Gammaj(1,:)==1),'.k','MarkerSize',ms,'LineWidth',lw);hold on;
    plot(index_time(Gammaj(2,:)==1),ksi_fem(Gammaj(2,:)==1),'.r','MarkerSize',ms,'LineWidth',lw);
    xlim(minmax(index_time))
    ylim([-0.6 0.6])
    datetick('x','yy','keeplimits');
    xlabel('Year','FontSize',fsize,'FontWeight','demi','Interpreter','LaTex')
    title(['Shape Parameter $\xi(U_t)$ for location ',infHandler.locNames{index_loc}],'Interpreter','LaTex','FontSize',fsize,'FontWeight','demi')
    
    h = legend('zero line','$\xi_1(U_t)$','$\xi_2(U_t)$', 'Orientation','horizontal');
    set(h,'Interpreter','LaTex','Box','on','FontSize',fsize,'FontWeight','demi')
    set(gca,'YTick',[-0.5 0 0.5],'YTickLabel',[-0.5 0 0.5],'FontSize',fsize);
    index_j = index_j + 1; 
end

end

