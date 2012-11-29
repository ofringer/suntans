%
% GET_COMPONENTS
%   [UAMP,UPHASE,VAMP,VPHASE,HAMP,HPHASE]=GET_COMPONENTS(YEAR,OMEGAT,LUN_NODE,V0U,AG)
%   calculates the tidal amplitudes and phases from the interpolated OTIS
%   data in the AG matrix.
%
%   This code has been adapted from Brian Dushaw's matlab scripts
%   obtained from http://909ers.apl.washington.edu/~dushaw/tidegui/tidegui.html
%
function [uamp,uphase,vamp,vphase,hamp,hphase]=get_components(YEAR,omegat,lun_node,v0u,AG)

  oneday=[335.62  0  322.55  1.97   334.63  0.99  -0.99  321.57];

  if YEAR < 1970 | YEAR > 2037,
    disp('constants for prediction year are not available')
    return;
  end

  years=lun_node(:,1);
  I=YEAR==years;
  vou=v0u(I,2:9);
  lunnod=lun_node(I,2:9);

  JJ=1:8;
  vou=(pi/180)*vou;
  oneday=(pi/180)*oneday;

  % Start with calculation of zonal currents.
  G=(pi/180)*AG(:,2)';
  A=AG(:,1)'.*lunnod;

  omegas = omegat(JJ)/24/3600;
  uamp = A(JJ);
  uphase = - oneday(JJ) + vou(JJ) - G;

  % Now get meridional currents.
  G=(pi/180)*AG(:,4)';
  A=AG(:,3)'.*lunnod;

  vamp = A(JJ);
  vphase = - oneday(JJ) + vou(JJ) - G;

  % Now get elevation.
  G=(pi/180)*AG(:,6)';
  A=AG(:,5)'.*lunnod;

  hamp = A(JJ);
  hphase = - oneday(JJ) + vou(JJ) - G;
