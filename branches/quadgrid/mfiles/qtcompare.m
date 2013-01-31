%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% qtcompare.m
% usage: plot interpolation results for tri and quad grid 
%        using lineplot function
% input: (just for lineplot) 
% 1) type = one of:
% i) 'q' nonhydrostatic pressure
% ii) 's' salinity
% iii) 'u','v','w' velocity fields
% iv) 'h free surface
% 2) datadir=location of data, i.e. '/home/data'
% 3) nstep= time step to plot
% 4) kevel=vertical level to plot
% 5) procs: total number of procs 
% 6) x,y: the line you want interpolate to 
% 7) METHOD: method for interpolation      
%    'nearest'   - Nearest neighbor interpolation
%    'linear'    - Linear interpolation (default)
%    'natural'   - Natural neighbor interpolation
%    'cubic'     - Cubic interpolation (2D only)
%    'v4'        - MATLAB 4 griddata method (2D only)
% 8) LINECOLOR: the color of the line
% 9) DIMENSION: 2 FOR 2D VECTOR(for simple line or y(x) is constant)
%               3 FOR 3D VECTOR
% 10) X: only for DIMENSION=2, 0 for x, 1 for y
% 11) the number of calculation steps
% output: comparsion figure and every figure for each step for making gif
% Made by Yun Zhang
% 01/25/2013 @Stanford
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nsteps=81; % the total number of calculation steps
gif=0;      % whether you want to save each plot for each step
movie=0;    % whether you want to plot all the plots for each step or  
            % just the step you want
PLOT='h';
n=21;        % the step you want to show
%dragorno=1;  % 1 with drag 0 no drag
METHOD='linear';
% directory of data
for dragorno=1:-1:0
if dragorno==1
    tridatadir='tridata/resultbc1-2/tridrag1amp0.01';
    trigriddir='tridata/resultbc1-2/tridrag1amp0.01';
    quaddatadir='quaddata/resultbc1-2/quaddrag1amp0.01';
    quadgriddir='quaddata/resultbc1-2/quaddrag1amp0.01';
    COLOR1='r';
    COLOR2='k';
else 
    tridatadir='tridata/resultbc1-2/trinodragamp0.01';
    trigriddir='tridata/resultbc1-2/trinodragamp0.01';
    quaddatadir='quaddata/resultbc1-2/quadnodragamp0.01';
    quadgriddir='quaddata/resultbc1-2/quadnodragamp0.01';
    COLOR1='r--';
    COLOR2='k--';
end
% assign the line
x=3:6:120;
y=30*ones(size(x));

if movie==0
   a11=lineplot(PLOT,trigriddir,tridatadir,n,1,0,x,y,METHOD,COLOR1,2,0); 
   hold on
   a22=lineplot(PLOT,quadgriddir,quaddatadir,n,1,0,x,y,METHOD,COLOR2,2,0);
   max(a11)/max(a22)
   legend('tri','quad')
end
end

% analytical solution part
T=40;
H=3.67;
k=pi/120;
amp=0.0367
if PLOT=='h'
    h1=amp*cos(k*x)*cos(2*pi/T*T/2);
    plot(x,h1,'g*-')
end

if PLOT=='w'
    h2=amp*cos(k*x)*cos(2*pi/T*T/4);
    w=-amp*2*pi/T/sinh(k*H)*(cosh(k*(h2+H))-cosh(k*(H-H)))./(h2+H).*cos(k*x)*sin(2*pi/T*T/4);
    plot(x,w,'g*-')
end

if PLOT=='u'
    h2=amp*cos(k*x)*cos(2*pi/T*T/4);
    u=amp*2*pi/T/(k*H)*sin(k*x)*sin(2*pi/T*T/4);
    plot(x,u,'g*-');
end
legend('tri with drag', 'quad with drag', 'tri no drag', 'quad no drag','analytical solution no drag')


if movie==1
    for i=1:nsteps
        lineplot(PLOT,tridatadir,trigriddir,i,1,0,x,y,METHOD,'r',2,0);
        hold on
        lineplot(PLOT,quaddatadir,quadgriddir,i,1,0,x,y,METHOD,'k',2,0); 
        legend('tri','quad')
        disp('press enter to continue!')
        hold off
        if gif==1
            a=num2str(i);
            saveas(gcf,a,'jpg');
        end
        pause
    end
end
        























