% FINDCLOSEST Find the processor and index of the grid cell closest
% to a given point.
%    [INDCLOSE,PROCCLOSE]=FINDCLOSEST(XC,YC,'cells.dat',32);
%    will return the index and the processor number of the Voronoi
%    point that is closest to the point XC, YC.
%
% Oliver Fringer
% Stanford University
% 12/20/06
function [indclose,procclose] = findclosest(xc,yc,fnameprefix,numprocs)
  
  mindist2 = 1e20;
  indclose = 0;
  procclose = 0;
  for proc=0:numprocs-1
    fname=sprintf('%s.%d',fnameprefix,proc);
    cells = load(fname,'-ascii');
    xv = cells(:,1);
    yv = cells(:,2);
    R2 = (xv-xc).^2 + (yv-yc).^2;
    dist2 = min(R2);
    if(dist2<mindist2)
      indclose = find(R2==dist2);
      procclose = proc;
      mindist2 = dist2;
    end
  end
    