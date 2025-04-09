figure('position',[100 500 700 500])
h = pcolor(x*1e6,y*1e6,abs(n));
xlabel('x (\mum)')
ylabel('y (\mum)')
colorbar
axis(1e6*[min(x) max(x) min(y) max(y)])
colormap(flipud(bone))
clim([min(min(n)) max(max(n))*1.5])