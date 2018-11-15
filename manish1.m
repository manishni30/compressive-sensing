
%length of the signal
N=1024;

%Number of random observations to take
K=256;

%Discrete frequency of two sinusoids in the input signal
k=20;

n=0:N-1;

%Sparse signal in frequency domain.
x=sin(2*pi*(k/N)*n); 
figure;
subplot(2,1,1);
plot(x)
grid on;
xlabel('Samples');
ylabel('Amplitude');
title('Original Signal,1024  samples with one frequency sinsuoids');

xf=fft(x);

xfmag=10*log10(abs(xf));

subplot(2,1,2);
plot(abs(xf))
grid on;
xlabel('Samples');
ylabel('Amplitude');
title('Frequency domain, 1024 coefficients with 2-non zero coefficients');

%creating dft matrix
B=dftmtx(N);
Binv=inv(B);

%Taking DFT of the signal
xf=B*x';

%Selecting random rows of the DFT matrix
q=randperm(N);

%creating measurement matrix
A=Binv(q(1:K),:);

%taking random time measurements
y=(A*xf);

%Calculating Initial guess
x0=A'*y;

%Running the recovery Algorithm
tic
xp=l1eq_pd(x0,A,[],y,1e-5);
toc

%recovered signal in time domain
xprec=real(Binv*xp);

figure;
subplot(2,1,1)
plot(abs(xf))
grid on;
xlabel('Samples');
ylabel('Amplitude');
title('Original Signal, Discrete Fourier Transform');

subplot(2,1,2)
plot(abs(xp),'r')
grid on;
xlabel('Samples');
ylabel('Amplitude');
title(sprintf('Recovered Signal, Discrete Fourier Transform sampled with %d samples',K));

figure;
subplot(2,1,1);
plot(x)
grid on;
xlabel('Samples');
ylabel('Amplitude');
title('Original Signal,1024 samples with one frequency sinsuoids');

subplot(2,1,2)
plot(xprec,'r')
grid on;
xlabel('Samples');
ylabel('Amplitude');
title(sprintf('Recovered Signal in Time Domain'));