%% signal initialization

n = 1/40000:1/40000:1/8;
f = sin(1394*pi*n) + sin(3266*pi*n);

%% Random Sampling

% b = Phi * f    (random sampling)
% c = Psi * f    where Psi = IDCT (sparsifying)

m = floor(rand(1,500)*length(f));
b = f(m);       % random samples of (f)
c = idct(f);    % sparsed (f)
plot(f,'b');hold
plot(m,b,'k.');title('Original Signal(f) and random sampless(b)');
axis([1,1000,-3,3])
legend('Original Signal','Randomly Sampled Signal')
figure, plot(c), axis([0,650,-10,10]);title('Sparse samples (c = IDCT(f))');

% f = Psi * c
% Psi = DCT
% f = DCT(c)
% c = IDCT(f)
%% solution 1 (x)

% Ax = b
% x = A\b

D = dct(eye(length(n),length(n)));
A = D(m,:);
% sound(f)
x = (A\b')';
b_hat = dct(x);
% sound(b_hat)

%% solution 2 (y)

% Ay = b
% y = pinv(A) * b

y = (pinv(A)*b')';
figure,plot(y),axis([0,650,-10,10])
figure,plot(dct(y)),axis([1,1000,-1.5,1.5])
%% solution 3 (s1)

s1 = l1eq_pd(y',A,A',b',5e-3,20); % L1-magic toolbox
