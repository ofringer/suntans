%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File name: plotsliceall.m
% Description: Extract data for a surface plot of all processors
%              from 2d suntans data.  This is an extension of plotslice.m
%              which only extracted data for one processor.  This
%              mfile extracts and merges the data for all processors.
%
% [x,z,data]=plotsliceall(plottype,dirname,timestep,numprocs);
% 
% plottype: 
%    'q': nonhydrostatic pressure
%    's': salinity
%    'u': u-velocity
%    'w': w-velocity
%    's0': background salinity
%    'h': free surface
%    'nut': eddy-viscosity
%    'kappat': scalar-diffusivity
%
%
% Oliver Fringer
% Stanford University
% 22 Jul 10
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,z,phi] = plotsliceall(type,datadir,step,nprocs)

  x = [];
  z = [];
  phi = [];

  for n=1:nprocs
    [x0,z0,phi0] = plotslice(type,datadir,step,n-1);
    x = [x,x0];
    z = [z,z0];
    phi = [phi,phi0];
  end

  x1 = x(1,:);
  [xs,is]=sort(x1);
  x = x(:,is);
  z = z(:,is);
  phi = phi(:,is);

  xd = diff(xs);
  ind = find(xd==0);
  x(:,ind) = [];
  z(:,ind) = [];
  phi(:,ind) = [];