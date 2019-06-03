
S_Rcoord = load('../R_swiss_coordinates.mat');
S_crd = load('../coordinates');

lon_km = (cat(1,S_crd.coordinates{2:N+1,4}));
lat_km = (cat(1,S_crd.coordinates{2:N+1,5}));
alt_km = (cat(1,S_crd.coordinates{2:N+1,6}));

scrsz = get(0,'ScreenSize');
figure1 = figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)],'Colormap',[0 0 0.5625;0 0 0.78125;0 0 1;0 0.0476190485060215 1;0 0.095238097012043 1;0 0.142857149243355 1;0 0.190476194024086 1;0 0.238095238804817 1;0 0.28571429848671 1;0 0.333333343267441 1;0 0.380952388048172 1;0 0.428571432828903 1;0 0.476190477609634 1;0 0.523809552192688 1;0 0.571428596973419 1;0 0.61904764175415 1;0 0.666666686534882 1;0 0.714285731315613 1;0 0.761904776096344 1;0 0.809523820877075 1;0 0.857142865657806 1;0 0.904761910438538 1;0 0.952380955219269 1;0 1 1;0.0625 1 0.9375;0.125 1 0.875;0.1875 1 0.8125;0.25 1 0.75;0.3125 1 0.6875;0.375 1 0.625;0.4375 1 0.5625;0.5 1 0.5;0.5625 1 0.4375;0.625 1 0.375;0.6875 1 0.3125;0.75 1 0.25;0.8125 1 0.1875;0.875 1 0.125;0.9375 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0;0.5 0 0]);
scnsize = get(0,'ScreenSize');

% plot swiss map
hold on;
h1 = surf(S_Rcoord.swiss_lon, S_Rcoord.swiss_lat, S_Rcoord.swiss_alt','EdgeColor','none','FaceLighting','phong');
shading interp  

% make city true
hold on;
plot3(lon_km/1000, lat_km/1000, 10*max(alt_km)*ones(1,N), '.k', 'MarkerSize', 40);
for j=1:17
    text(lon_km(j)/1000+8, lat_km(j)/1000, 10*max(alt_km), S_crd.coordinates{j+1,1},'HorizontalAlignment','left','FontSize',26,'FontWeight','demi')
end