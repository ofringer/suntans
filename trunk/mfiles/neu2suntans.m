% NEU2SUNTANS  Convert Gambit .neu file to suntans readable grid.
%   NEU2SUNTANS(MESH_FILE,DATADIR) reads the .neu file specified
%   by MESH_FILE and creates the points.dat, edges.dat, cells.dat files
%   that are used by SUNTANS in the directory specified by DATADIR.
%   If only one argument is specified then it is assumed that DATADIR is
%   the local directory.  If no arguments are specified, then the
%   default neu file is assumed to be default.neu

% Sheng Chen
% ICME, Stanford University
% Febrary 02, 2006
function neu2suntans(varargin)
  
  if(nargin==0)
    mesh_file='default.neu';
    datadir='.';
  elseif(nargin==1)
    mesh_file=varargin{1};
    datadir='.';
  elseif(nargin==2)
    mesh_file=varargin{1};
    datadir=varargin{2};
  else
    error('neu2suntans must be called with one or two arguments.');
  end

  if(~exist(mesh_file,'file'))
    error(sprintf('File %s does not exist',mesh_file));
  end
  if(~exist(datadir,'dir'))
    error(sprintf('Directory %s does not exist',datadir));
  end
  points_file = [datadir,'/points.dat'];
  edges_file = [datadir,'/edges.dat'];
  cells_file = [datadir,'/cells.dat'];

  % Start the counter
  tic;

  % read data from Gambit Neutral file
  % you can change the input file name if needed
  disp(sprintf('Processing %s ...', mesh_file))
  meshf = fopen(mesh_file,'r');
  
  info = fscanf(meshf,'%s',22);
  num = fscanf(meshf,'%d'); np = num(1); nc = num(2);
  
  info = fscanf(meshf,'%s',4);
  point = fscanf(meshf,'%d %e %e',[3, np]); point = point';
  
  info = fscanf(meshf,'%s',3);
  cell = fscanf(meshf,'%d %d %d %d %d %d',[6, nc]); cell = cell';
  
  % sort the point number such that in the adjacency matrix
  % we only need to care about the upper triangular part
  cell = sort(cell(:, 4:6), 2);
  
  info = fscanf(meshf,'%s',13);
  group = fscanf(meshf,'%d');
  
  info = fscanf(meshf,'%s',6);
  num1 = fscanf(meshf,'%d',4);
  if(isempty(num1))
    nBC1 = 0;
    BC1 = [];
  else
    nBC1 = num1(2);
    BC1 = fscanf(meshf,'%d %d %d',[3, nBC1]); BC1 = BC1';
  end
  
  info = fscanf(meshf,'%s',6);
  num2 = fscanf(meshf,'%d',4);
  if(isempty(num2))
    nBC2 = 0;
    BC2 = [];
  else
    nBC2 = num2(2);
    BC2 = fscanf(meshf,'%d %d %d',[3, nBC2]); BC2 = BC2';
  end
  
  status = fclose(meshf);
  
  % convert Gambit to Suntans
  % x0 and y0 are the values to be subtracted in order to keep the origin at (0,0)
  x = point(:, 2); y = point(:, 3);
  points = [x, y, zeros(np, 1)];
  
  % adjacency matrix
  A = sparse(np, np);
  Vp1 = sparse(np, np); Vp2 = sparse(np, np);
  xv = zeros(nc, 1); yv = zeros(nc, 1);

  for i = 1:nc
    % calculate the circumcenter of each cell, which is the Voronoi point in
    % the triangulation
    
    p1 = cell(i,1); p2 = cell(i,2); p3 = cell(i,3);
    [xv(i), yv(i)] = circumcenter(x(p1), y(p1), x(p2), y(p2), x(p3), y(p3));

    if A(p1, p2) == 0
      Vp1(p1, p2) = i;
    else 
      Vp2(p1, p2) = i;
    end %if
    
    if A(p1, p3) == 0
      Vp1(p1, p3) = i;
    else 
      Vp2(p1, p3) = i;
    end %if
    
    if A(p2, p3) == 0
      Vp1(p2, p3) = i;
    else 
      Vp2(p2, p3) = i;
    end %if
    
    A(p1, p2) = A(p1, p2) + 1; A(p1, p3) = A(p1, p3) + 1; A(p2, p3) = A(p2, p3) + 1;
  end
  
  [I,J] = find(A>0);
  ne = length(I);
  % edge = full([I-1, J-1, A((J-1)*np+I)~=2, Vp1((J-1)*np+I)-1, Vp2((J-1)*np+I)-1]);
  % using the following way to avoiding the index over the integer limit
  edge = zeros(ne, 5);
  for k = 1:ne
    i = I(k); j = J(k); edge(k, :) = [i-1, j-1, A(i,j)~=2, Vp1(i,j)-1, Vp2(i,j)-1];
  end %k
  
  % apply offshore BC and set the corresponding type to 2
  for i = 1:nBC2
    for j = 1:ne
      if edge(j, 5) == -1 & edge(j, 4) == BC2(i,1) - 1
        edge(j, 3) = 2;
      end    
    end
  end    
  
  % obtain the neighbor information of each cell
  neighbor = zeros(nc,3);
  for i = 1:nc
    p1 = cell(i,1); p2 = cell(i,2); p3 = cell(i,3);
    if Vp1(p1,p2) == i
      neighbor(i, 1) = Vp2(p1,p2);
    else
      neighbor(i, 1) = Vp1(p1,p2);
    end
    
    if Vp1(p1,p3) == i
      neighbor(i, 2) = Vp2(p1,p3);
    else
      neighbor(i, 2) = Vp1(p1,p3);
    end
    
    if Vp1(p2,p3) == i
      neighbor(i, 3) = Vp2(p2,p3);
    else
      neighbor(i, 3) = Vp1(p2,p3);
    end
  end
  
  % write tada to Suntans required files
  pointsf = fopen(points_file, 'w');
  fprintf(pointsf, '%12.10e %12.10e %d\n', points');
  status = fclose(pointsf);
  
  edgesf = fopen(edges_file, 'w');
  fprintf(edgesf, '%1.0f %1.0f %1.0f %1.0f %1.0f\n', edge');
  status = fclose(edgesf);
  
  cells = [xv, yv, cell-1, neighbor-1];
  cellsf = fopen(cells_file,'w');
  fprintf(cellsf, '%12.10e %12.10e %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f\n', cells');
  status = fclose(cellsf);        
  
  fprintf('Execution time = %.3f min (%.1f sec)\n',toc/60,toc);
  