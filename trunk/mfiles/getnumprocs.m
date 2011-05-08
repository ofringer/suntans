% 
% GETNUMPROCS Determine number of processors in a suntans run.
%   NP = GETNUMPROCS(datadir) returns the number of processors in a
%   suntans run with data in datadir.  Returns -1 if no data exists
%   in the directory or if the directory does not exist.
% 
function np = getnumprocs(datadir)

  suntansfile=[datadir,'/suntans.dat'];
  if(~exist(suntansfile,'file'))
    error(sprintf('%s does not exist',suntansfile));
    np=-1;
    return;
  end
  
  fvals = {'FreeSurfaceFile',...
           'HorizontalVelocityFile',...
           'VerticalVelocityFile',...
           'SalinityFile',...
           'BGSalinityFile',...
           'TemperatureFile',...
           'PressureFile'};

  for fval=fvals
    root=getvalue(suntansfile,fval);
    np=0;
    while(true)
      if(~exist(sprintf('%s/%s.%d',datadir,root,np)))
        break;
      else
        np=np+1;
      end
    end
  end

  if(np==0)
    np=-1;
  end