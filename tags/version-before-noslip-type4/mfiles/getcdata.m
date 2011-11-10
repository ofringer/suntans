%
% GETCDATA
% DATA = GETCDATA(FID,ARRAYSIZE,STEP,PRECISION);
%   extracts binary data from the specified file descriptor.  As
%   an example, if a file contains double precision arrays of length
%   N=10, the n=10th array is extracted with:
%   
%   data = getcdata(fid,N,n,'float64');
%
% Oliver Fringer
% Stanford University
% 27 July 05
%
function data = getcdata(fid, arraysize, step, precision)

  tic;
  if(fid>0),
    if(precision=='float32'),
      space = 4*arraysize;                    
    elseif(precision=='float64')
      space = 8*arraysize;
    end

    frewind(fid);
    fseek(fid,(step-1)*space,0);
    data = fread(fid,arraysize,precision);
  else
    fprintf('Error in function: getdata...undefined fid\n');
    data = zeros(arraysize,1);
  end
  
  fprintf('Read data at a rate of %.2f Mb/sec\n',space/1e6/toc);
