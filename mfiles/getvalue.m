function value = getvalue(fname,str)
  
  fid = fopen(fname,'r');
  if(fid==-1)
    fprintf('Error opening %s\n',fname);
    return;
  end
  
  while(1)
    line = fgetl(fid);
    if(line==-1)
      break;
    else
      if(line(1)~='#')
        strs=strread(line,'%s');
        if(strcmp(strs{1},str))
          value=str2num(strs{2});
          break;
        end
      end
    end
  end
  
  