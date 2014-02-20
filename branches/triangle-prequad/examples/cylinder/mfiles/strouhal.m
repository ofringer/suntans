%
% Calculate the Strouhal number from profile data obtained with the
% cylinder example.  This run gives a Strouhal number of 0.25.
%
addpath '../../../../mfiles';
dt0 = str2num(getvalue('../data/suntans.dat','dt'));
ntoutProfs = str2num(getvalue('../data/suntans.dat','ntoutProfs'));

dt = ntoutProfs*dt0;
d = 0.1;
u0 = 1;

datadir='../data';
fname = [datadir,'/u.dat.prof'];

fid = fopen(fname,'rb');
udata = fread(fid,'float64');
u = udata(1:3:end);
v = udata(2:3:end);
w = udata(3:3:end);

N = length(v);
t = dt*[1:N];
omega = 2*pi*[-N/2:N/2-1]/t(end);

Vhat = fftshift(fft(v))/N;

Vamp = sqrt(Vhat.*conj(Vhat));
Vphase = atan(imag(Vhat)./real(Vhat));

peak = find(Vamp==max(Vamp));

Vpeak = Vamp(peak);
omegapeak = omega(peak);
Vim = Vhat(peak);

omegapeak = omegapeak(end);
Vpeak = Vpeak(end);
Vim = Vim(end);

figure(1);
subplot(2,1,1)
loglog(omega(N/2+1:end)/2/pi,Vamp(N/2+1:end)/u0,'k-',...
       omegapeak/2/pi,Vpeak,'ko');
text(omegapeak/2/pi+1,Vpeak,sprintf('Strouhal number is f d/U = %.2f\n',(omegapeak/2/pi)*d/u0));

xlabel('Frequency (Hz)');
ylabel('v/u_0');

subplot(2,1,2)
plot(omegapeak*t/2/pi,v/u0,'k-')
xlabel('f*t');
ylabel('v/u_0');




