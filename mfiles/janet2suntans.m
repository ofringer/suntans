%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% janet2suntans.m for hybrid version suntans
% usage:
% 1. add 1 columns for nfaces
% 2. find the missing point and face for cells.dat
% 3. output the complete files for suntans
% 4. can work with any hybrid grid from janet
% 5. change edge type
% input:
% 1. any suntans type data to add nfaces column
% 2. Janet hybrid gird (nontri & tri)with edge depth and
% output in suntans type
% shortcomings: maybe not correct in order of neigh for
% cells with nfaces>4
% Made by Yun Zhang
% 01/22/2013 @Stanford
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read original data from points cells edges
pointdata = load(['points.dat']);
celldata = load(['cells.dat']);
edgedata = load(['edges.dat']);
% change edgedata or just create new celldata only?
edgechange=0; % 0 for celldata only 

% x y for points
px=pointdata(:,1);
py=pointdata(:,2);
% x y for cell voronoi points
cx=celldata(:,1);
cy=celldata(:,2);
% marker for boundary
mark=edgedata(:,3);
nc=length(cx);
ne=length(mark);
nfaces=zeros(nc,1);
% index for edge points, change to 1-based from zero-based
pi=edgedata(:,1:2);
% index for edges, change to 1-based from zero-based 
% (don't change the number for boundary edges)
ci=edgedata(:,4:5);

% here we can change the marker and create new edge
if edgechange==1
    mark(find(mark==2))=1; % add whatever you want
end

% start from edge
for i=1:ne
    for j=1:2
        if ci(i,j)>=0
            k=ci(i,j)+1;
            nfaces(k)=nfaces(k)+1;
        end
    end
end

% locate where non-triangular grids are
loc=find(nfaces~=3);

% % find the different number of edges except 3
% numedge=unique(nfaces);
% unqedge=numedge(find(numedge~=3));
% unqnum=length(unqedge);

% maxfaces: the biggest number of edges for a single cell
maxfaces=max(nfaces);

% points index which make up of cells
cp=-10*ones(nc,maxfaces);
cp(:,1:3)=celldata(:,3:5);

% neighbouring cell 
cneigh=-10*ones(nc,maxfaces);
cneigh(:,1:3)=celldata(:,6:8);

order=[2,1];
for i=1:length(loc)
    kp=4;
    kn=4;
    kbc=0;
    for j=1:2
        % find edges which make nontriangular grids
        loc2=find(ci(:,j)==(loc(i)-1));
        % add neigh and point
        for k=1:length(loc2)
            % add point information
            ep=pi(loc2(k),:);
            check1=ismember(ep,cp(loc(i),:));
            check2=find(check1==0);
            if length(check2)~=0
                cp(loc(i),kp:(kp+length(check2)-1))=ep(check2);
                kp=kp+length(check2);
            end
            
            % add neigh information
            ec=ci(loc2(k),order(j));
            if ec==-1
                kbc=kbc+1;
            end    
            check1=ismember(ec,cneigh(loc(i),1:3));
            if check1==0
                cneigh(loc(i),kn)=ec;
                kn=kn+1;
            end   
        end
    end
    loc3=find(cneigh(loc(i),:)==-1);
    if length(loc3)~=kbc
        % maybe not right in order of neigh when nfaces>5 at boundary cell
        % in which more than 2 edges are boundary edges (neigh=-1)
        for kkk=1:(kbc-length(loc3))
            cneigh(loc(i),kn)=-1;
            kn=kn+1;
        end
    end
end

% check results
% show edges
for i=1:ne
    plot(px(pi(i,:)+1),py(pi(i,:)+1));
    hold on
end

% show points and neigh for cells you assign
npp=input('which cell you want to see in the whole grid?  input -1 to exit               ');
while npp~=-1
if npp>nc
    disp('beyong the limit')
else
    ncc=cneigh(npp+1,1:nfaces(npp+1))+1;
    ncc(find(ncc==0))=[];
    plot(cx(ncc),cy(ncc),'r*','linewidth',1.5); 
    plot(px(cp(npp+1,1:nfaces(npp+1))+1),py(cp(npp+1,1:nfaces(npp+1))+1),'k*','linewidth',1.5);
end
npp=input('which cell you want to see in the whole grid?  input -1 to exit               ');
end



% output new cell data
 cells_file = ['cells_new.dat'];
 cellsf = fopen(cells_file,'w');
 for i=1:nc
     fprintf(cellsf, '%1.0f %12.10e %12.10e ', nfaces(i),cx(i),cy(i));
     for j=1:nfaces(i)
         fprintf(cellsf, '%1.0f ', cp(i,j));
     end
     for j=1:(nfaces(i)-1)
         fprintf(cellsf, '%1.0f ', cneigh(i,j));
     end
     fprintf(cellsf, '%1.0f\n',cneigh(i,nfaces(i)));
 end  
 status = fclose(cellsf);       

% output new edge data
if edgechange==0
edges_file=['edges_new.dat'];
edgesf = fopen(edges_file,'w');
for i=1:ne
    fprintf(edgesf, '%1.0f %1.0f %1.0f %1.0f %1.0f\n', pi(i,1), pi(i,2), mark(i), ci(i,1), ci(i,2));
end
status = fclose(edgesf)
end

















