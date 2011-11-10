function alldata = readalldata(dirname,fname,component,type,numprocs,nstep,Nkmax,Nc,i,repeats,VERBOSE)

  EMPTY=999999;
  
  if(type=='u')
    numcomponents=3;
    if(component<1 | component>3)
      fprintf('Error in function readalladata: \n');
      fprintf('\t component must be 1 (for u), 2 (for v), or 3 (for w).\n');
      alldata = -1;
      return;
    end
  elseif(type=='s')
    numcomponents=1;
    component=1;
  else
    fprintf('Error in function readalladata: \n');
    fprintf('\t type variable must be either ''u'' or ''s''.\n');
    alldata=-1;
    return;
  end

  nn=1;
  if(VERBOSE)
    fprintf('Loading data in %s.* files...\n',fname);
  end
  for n=0:numprocs-1
    datafile = [dirname,'/',fname,'.',num2str(n)];
    datafid = fopen(datafile,'rb');
    if(datafid<0)
      fprintf('Error opening %s.\n',datafile);
      alldata=-1;
      return;
    end
    arraysize=Nc(nn)*Nkmax*numcomponents;
    fseek(datafid,(nstep-1)*8*arraysize,0);
    data = fread(datafid,arraysize,'float64');
    if(length(data)~=numcomponents*Nc(nn)*Nkmax)
      fprintf('Error reading %s (maybe needs to be a velocity data file?)\n',datafile);
      alldata=-1;
      return;
    end
    if(numcomponents==1)
      data = reshape(data,Nc(nn),Nkmax);
      datapp{nn} = data;
    else
      data = reshape(data,Nc(nn),numcomponents,Nkmax);      
      datapp{nn} = squeeze(data(:,component,:));
    end
    
    nn=nn+1;
  end
  
  if(VERBOSE)
    fprintf('Merging processors...\n');
  end
  for n=0:numprocs-1
    if(n==0)
      alldata = datapp{1};
    else
      alldata(end+1:end+Nc(n+1),:) = datapp{n+1};
    end
  end
  
  if(VERBOSE)
    fprintf('Removing ghost points...\n');
  end
  alldata(i(repeats),:) = [];
  alldata(find(alldata==EMPTY))=nan;
  