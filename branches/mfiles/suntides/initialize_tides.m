%
% INITIALIZE_TIDES 
%   INITIALIZE_TIDES(PATH) Load OTIS tidal data from the given path
%   and makes the variables global for use by the get_tides function.
%
%   Note that only allconsts.mat and latlons.mat are required to exist
%   in the directory defined by PATH.
%
%   This code has been adapted from Brian Dushaw's matlab scripts
%   obtained from http://909ers.apl.washington.edu/~dushaw/tidegui/tidegui.html
%
function initialize_tides(tidespath)

  global lun_node v0u lat lon latu lonu latv lonv h ...
      ui1 ur1 vi1 vr1 eli1 elr1 ui2 ur2 vi2 vr2 eli2 elr2 ...
      ui3 ur3 vi3 vr3 eli3 elr3 ui4 ur4 vi4 vr4 eli4 elr4 ...
      ui5 ur5 vi5 vr5 eli5 elr5 ui6 ur6 vi6 vr6 eli6 elr6 ...
      ui7 ur7 vi7 vr7 eli7 elr7 ui8 ur8 vi8 vr8 eli8 elr8
  
  %load in all the tide data and the model bathymetry - 79MB
  load([tidespath,'/allconsts.mat']);
  % now have h, elr, eli, ur, ui, vr, vi defined for all 8 constituents.
  
  %load in the lat and lon values for h, el, u, v
  load([tidespath,'/latlons.mat']);
  %now have lat,lon; latu, lonu; latv; lonv

  str = which('suntides.m');
  datapath = str(1:findstr(str,'suntides.m')-2);

  load([datapath,'/data/lun_node.dat']);
  load([datapath,'/data/v0u.dat']);
  
