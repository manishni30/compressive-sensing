function x1 = resconstruct(y,n,cluster_id) 

    %___MEASUREMENT MATRIX___
    m = length(y); % NOTE: small error still present after increasing m to 1500;
    
        %__ALTERNATIVES TO THE ABOVE MEASUREMENT MATRIX___ 
        %Phi = (sign(randn(m,n))+ones(m,n))/2; % micro mirror array (mma) e.g. single
        %pixel camera Phi = orth(Phi')'; % NOTE: See Candes & Romberg, l1
        %magic, Caltech, 2005.

    %___COMPRESSION___
    rng(cluster_id);
    Phi = randn(m,n);

    %___THETA___
    % NOTE: Avoid calculating Psi (nxn) directly to avoid memory issues.
    Theta = zeros(m,n);
    for ii = 1:n
   
        ek = zeros(1,n);
        ek(ii) = 1;
        psi = idct(ek)';
        Theta(:,ii) = Phi*psi;
    end
    
    %___l2 NORM SOLUTION___ s2 = Theta\y; %s2 = pinv(Theta)*y
    s2 = pinv(Theta)*y;


    %___BP SOLUTION___
    s1 = l1eq_pd(s2,Theta,Theta',y,5e-3,20); % L1-magic toolbox
    %x = l1eq_pd(y,A,A',b,5e-3,32);

   
    %___IMAGE RECONSTRUCTIONS___
    x2 = zeros(n,1);
    for ii = 1:n
   
        ek = zeros(1,n);
        ek(ii) = 1;
        psi = idct(ek)';
        x2 = x2+psi*s2(ii);
    end

    x1 = zeros(n,1);
    
    for ii = 1:n
    
        ek = zeros(1,n);
        ek(ii) = 1;
        psi = idct(ek)';
        x1 = x1+psi*s1(ii);
    end
    
end


