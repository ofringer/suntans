%
% PRINTMOVIE(15,'/tmp/movies/movie',3,50) will print a jpeg and a gif 
%    of the 15th frame of a series of images in figure 3 into the 
%    directory /tmp/movies, and the filename will be called 
%    '/tmp/movies/movie0015.jpg'.  The jpeg will be scaled by 50% 
%    to form the gif file '/tmp/movies/movie0015.gif'.  Numbering is
%    assigned so that movie2.gif --> movie0002.gif, which comes before
%    movie10.gif and the movie files are appropriately ordered for use with
%    gifsicle or convert.
%
% Oliver Fringer
% Stanford University
% 10 July 2001
%
function printmovie(n,fname,cf,scale)
  
  figure(cf);
  
  if(n<10) 
    str = ['000',num2str(n)];
  elseif(n>9 & n <100)
    str = ['00',num2str(n)];
  elseif(n>99 & n<1000)
    str = ['0',num2str(n)];
  else
    str = num2str(n);
  end
 
  eval(['print -depsc2 ',fname,str]);
  unix(['convert -crop 0x0 -geometry ',sprintf('%d%%',scale),' ',fname,str,'.eps ',fname,str,'.gif']);
