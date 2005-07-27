addpath('../../../mfiles');

datadir='../../../main/examples/iwaves/data';

[x,y,z]=plotslice(5,datadir,10,0);

clf;
pcolor(x/1000,y,z);
shading flat;
daspect([1 75 1]);
hc = colorbar('horiz');
pc = get(hc,'position');
set(hc,'position',[pc(1) pc(2)+1.75*pc(4) pc(3) pc(4)]);

xlabel('x (km)');
ylabel('y (m)');

