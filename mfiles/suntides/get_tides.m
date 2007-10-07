%
% GET_TIDES Retrieve tidal constituents from OTIS tidal data.
%   [omegas,uamp,uphase,vamp,vphase,hamp,hphase]=get_tides(YEAR,latw,lonw)
%   returns the top 8 tidal constituents given the year, latitude
%   (North), and longitude (East), in a way that the velocity (or
%   elevation) at that location and time (in seconds since the
%   beginning of the year) can be obtained with:
%
%   u = 0;
%   v = 0;
%   h = 0;
%   for n=1:length(omegas)
%     u = u + uamp(n)*cos(omegas(n)*time + uphase(n));
%     v = v + vamp(n)*cos(omegas(n)*time + vphase(n));
%     h = h + hamp(n)*cos(omegas(n)*time + hphase(n));
%   end
%
%   This code has been adapted from Brian Dushaw's matlab scripts
%   obtained from http://909ers.apl.washington.edu/~dushaw/tidegui/tidegui.html.
%

%
%   Oliver Fringer
%   Stanford University
%   9 Sep 07
%
function [omegas,uamp,uphase,vamp,vphase,hamp,hphase]=get_tides(YEAR,latw,lonw)

  global lun_node v0u lat lon latu lonu latv lonv h ...
      ui1 ur1 vi1 vr1 eli1 elr1 ui2 ur2 vi2 vr2 eli2 elr2 ...
      ui3 ur3 vi3 vr3 eli3 elr3 ui4 ur4 vi4 vr4 eli4 elr4 ...
      ui5 ur5 vi5 vr5 eli5 elr5 ui6 ur6 vi6 vr6 eli6 elr6 ...
      ui7 ur7 vi7 vr7 eli7 elr7 ui8 ur8 vi8 vr8 eli8 elr8
  
  ipred=1;
  
  AG=[];
  % Tidal periods in hours
  tm2=12.4206;
  ts2=12.;
  tn2=12.6583;
  tk2=11.9672;
  to1=25.8193;
  tk1=23.9345;
  tp1=24.0659;
  tq1=26.8724;

  % omegat has the tidal frequencies in rad/yearday.   
  % omegas has them in rad/s
  omegat=48.*pi./[tm2 ts2 tn2 tk2 to1 tk1 tp1 tq1];
  omegas=omegat/(3600*24);

  uamp = [];
  uphase = [];
  vamp = [];
  vphase = [];
  hamp = [];
  hphase = [];
  
  if(YEAR==0)
    return;
  end
  
  % the order of frequencies is M2, S2, N2, K2, O1, K1, P1, Q1
  % ui,ur,vi,vr,eli,elr are all loaded in predict.m
  for iop=1:8,
    if iop==1,
      ui=ui1;
      ur=ur1;
      vi=vi1;
      vr=vr1;
      eli=eli1;
      elr=elr1;
    elseif iop==2,
      ui=ui2;
      ur=ur2;
      vi=vi2;
      vr=vr2;
      eli=eli2;
      elr=elr2;
    elseif iop==3,
      ui=ui3;
      ur=ur3;
      vi=vi3;
      vr=vr3;
      eli=eli3;
      elr=elr3;
    elseif iop==4,
      ui=ui4;
      ur=ur4;
      vi=vi4;
      vr=vr4;
      eli=eli4;
      elr=elr4;
    elseif iop==5,
      ui=ui5;
      ur=ur5;
      vi=vi5;
      vr=vr5;
      eli=eli5;
      elr=elr5;
    elseif iop==6,
      ui=ui6;
      ur=ur6;
      vi=vi6;
      vr=vr6;
      eli=eli6;
      elr=elr6;
    elseif iop==7,
      ui=ui7;
      ur=ur7;
      vi=vi7;
      vr=vr7;
      eli=eli7;
      elr=elr7;
    elseif iop==8,
      ui=ui8;
      ur=ur8;
      vi=vi8;
      vr=vr8;
      eli=eli8;
      elr=elr8;
    end
    
    % WE WILL ONLY INTERPOLATE USING DATA WITHIN 2 GRID POINTS
    % note there is potentially a problem near land because land has
    % a value of zero for the tide - the interpolation will work using
    % this value of zero.  Therefore, results within a grid point of
    % land may be suspect.
    %
    DEL_LON=2*(lonu(1,2)-lonu(1,1));
    DEL_LAT=2*(latu(2,1)-latu(1,1));
    I1=((lonw-DEL_LON) < lonu(1,:) & lonu(1,:) < (lonw+DEL_LON)) ;
    I2=((latw-DEL_LAT) < latu(:,1) & latu(:,1) < (latw+DEL_LAT));

    % Interpolation method
    interp_method = '*linear';
    
    if sum(I1)>=3 & sum(I2)>=3,
      lon_tmp=lonu(I2,I1);
      lat_tmp=latu(I2,I1);
      ur_tmp=ur(I2,I1);
      ui_tmp=ui(I2,I1);
      UR=interp2(lon_tmp,lat_tmp,ur_tmp,lonw,latw,interp_method);
      UI=interp2(lon_tmp,lat_tmp,ui_tmp,lonw,latw,interp_method);
      Au14=sqrt(UR.^2 + UI.^2);
      Gu14=atan2(-UI,UR)*(180/pi);
      
      I1=((lonw-DEL_LON) < lonv(1,:) & lonv(1,:) < (lonw+DEL_LON)) ;
      I2=((latw-DEL_LON) < latv(:,1) & latv(:,1) < (latw+DEL_LON));
      lon_tmp=lonv(I2,I1);
      lat_tmp=latv(I2,I1);
      vr_tmp=vr(I2,I1);
      vi_tmp=vi(I2,I1); 
      VR=interp2(lon_tmp,lat_tmp,vr_tmp,lonw,latw,interp_method);
      VI=interp2(lon_tmp,lat_tmp,vi_tmp,lonw,latw,interp_method);
      Av14=sqrt(VR.^2 + VI.^2);
      Gv14=atan2(-VI,VR)*(180/pi); 
      
      I1=((lonw-DEL_LON) < lon(1,:) & lon(1,:) < (lonw+DEL_LON)) ;
      I2=((latw-DEL_LAT) < lat(:,1) & lat(:,1) < (latw+DEL_LAT));
      lon_tmp=lon(I2,I1);
      lat_tmp=lat(I2,I1);
      h_tmp=h(I2,I1);
      elr_tmp=elr(I2,I1);
      eli_tmp=eli(I2,I1);
      h14=interp2(lon_tmp,lat_tmp,h_tmp,lonw,latw,interp_method);
      Ar=interp2(lon_tmp,lat_tmp,elr_tmp,lonw,latw,interp_method);
      Ai=interp2(lon_tmp,lat_tmp,eli_tmp,lonw,latw,interp_method);
      
      A=100*sqrt(Ar.^2 + Ai.^2);
      % minus sign below because Gary's sign convention for phase is opposite
      % to the rest of the world's.
      G=(180/pi)*atan2(-Ai,Ar);
      
      % the values of Au14 and Av14 are actually transports in m^2/s.  Divide by
      % the model depth and multiply by 100 to convert to current in
      % (cm/s)
      Au14=100*Au14./h14;
      Av14=100*Av14./h14;
      Gu14( Gu14 > 360)=Gu14( Gu14 > 360)-360;
      Gv14( Gv14 > 360)=Gv14( Gv14 > 360)-360;
      Gu14( Gu14 <   0)=Gu14( Gu14 < 0  )+360;
      Gv14( Gv14 <   0)=Gv14( Gv14 < 0  )+360;
      G( G > 360)=G( G > 360)-360;
      G( G <   0)=G( G < 0  )+360;
      
      AG=[AG;[Au14 Gu14 Av14 Gv14 A G]];
    else
      ipred=0;
    end
  end
  
  if ipred==1,
    % now pass AG to get_components
    [uamp,uphase,vamp,vphase,hamp,hphase]=get_components(YEAR,omegat,lun_node,v0u,AG);
  else
    fprintf('\n')
    fprintf('ERROR:  Point %.2f E,%.2f N is out of the range of tidal data:\n',lonw,latw);
    fprintf('\t\t Range is: (%.2f to %.2f E, %.2f to %.2f N)\n',...
            min(min(lon)),max(max(lon)),min(min(lat)),max(max(lat)));
    fprintf('Make sure lon/lat coordinates are E/N!\n');
    fprintf('\n')
    return;
  end
  
  