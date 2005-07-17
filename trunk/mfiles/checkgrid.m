%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CHECKGRID Check Voronoi distances on an unstructured grid.
%   CHECKGRID(DGMIN,NBINS,EDGESFILE,CELLSFILE) checks an 
%   unstructured grid defined by the EDGESFILE and CELLSFILE
%   for Voronoi distances (distances between Voronoi points on either 
%   side of an edge) that are less than DGMIN, and displays those edges in
%   red on a plot of the Voronoi points.  This function also plots
%   a histogram of the Voronoi distances with the number of bins
%   as NBINS.
%  
%   Oliver Fringer
%   Stanford University
%   16 July 05
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkgrid(dgmin,nbins,edgesfile,cellsfile)

  e = load(edgesfile,'-ascii');
  c = load(cellsfile,'-ascii');

  xv = c(:,1);
  yv = c(:,2);

  g2j = 1+e(:,4);
  g2jp1 = 1+e(:,5);

  g2jp1(find(g2j==0))=[];
  g2j(find(g2j==0))=[];
  g2j(find(g2jp1==0))=[];
  g2jp1(find(g2jp1==0))=[];

  Dg = sqrt((xv(g2j)-xv(g2jp1)).^2+(yv(g2j)-yv(g2jp1)).^2);
  ind = find(Dg<dgmin);

  figure(1);
  hist(Dg,linspace(min(Dg),max(Dg),nbins));
  xlabel('Voronoi Distance');
  ylabel('Number of edges');
  title('Voronoi Distance Histogram');

  figure(2);
  clf;
  hold on;
  plot(xv,yv,'k.','markersize',1);
  plot([xv(g2j(ind)),xv(g2jp1(ind))]',...
       [yv(g2j(ind)),yv(g2jp1(ind))]','r.-','markersize',10);
  axis image;
  xlabel('x');
  ylabel('y');
  title(sprintf('Edges with length less than %.2f shown in red',dgmin));
