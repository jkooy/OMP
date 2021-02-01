clc;clear;

%%%%%%%%%%%%%%%%%%%%%%% (d) Noisy case: (n != 0) %%%%%%%%%%%%%%%%%%%%%%%%

% Generate different phase transition plots for the following values of N: 20, 50 and 100.
N_list=[20,50,100];

%%(a)
    sigma_list = [1e-3,1e-1]; 
    % two values of s (one small and one large)
    N=N_list(n);
    M_max=0.8*N;
    M_max=35;
    % smax is chosen to be a reasonably large integer which is smaller than N
    s_max=0.3*N;
    s_max=15;
   
    
    ESR_matrix = zeros(M_max,s_max);
    aNE_matrix = zeros(M_max,s_max);
    
    for M = 1:M_max
        for s = 1:s_max 
             
            ESR=0;
            aNE=0;
            for i=1:times
%% a
                
%%%%%%%%%%%%%%%%%%%%%%%(b)Experimental setup:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
                % Generate A as a random matrix with independent and identically distributed entries drawn from the standard normal distribution

                A = randn(M,N);

                % Normalize the columns of A, not use in-built library functions.
                A_vn = sqrt(diag(A'*A))';
                A = A./A_vn;         
                % norm(A(:,1))

%% b    
                % Generate the sparse vector x with random support of cardinality s

        %       s indices are
        %       generated uniformly at random from integers 1 to N
                S = randi([1 N],s,1);
                S = sort(S); % sort S
        %       non-zero entries drawn as uniform random variables
                e = zeros(s,1);
                for j = 1:s
                    U = randi([0 1]);
                    if U==1
                        e(j) = 1 + 9*rand;
                    elseif U==0
                        e(j) = -1 + (-9)*rand;
                    end
                end
                % generate x
                x = zeros(N,1);
                x(S) = e;
%% c 
                noise = randn(M,1);
                % Noiseless case
                y = A(:,S)*e;
                
                [x_hat,S_hat] = OMP1(A,y);
                              
                if length(S_hat) == s
                    if sum(sort(S_hat)==S)==s
                        ESR = ESR+1;
                    end
                end
           
                % identify ANE(average Normalized Error)
                aNE = aNE + PM(x,x_hat);                
            end
            ESR_matrix(M,s)=ESR/times;
            aNE_matrix(M,s)=aNE/times;
        end
        
        disp(M); 
    end
    figure(1);
    title('ESR');
    imagesc(ESR_matrix);
    colormap gray
 
    figure(2);
    title('aNE');
    imagesc(aNE_matrix);
    colormap gray
end




function [x_hat, S_hat] = OMP1(A,y)
N = size(A,2);
S = zeros(N,1);
r = y;
r0 = 100;
k = 1;
while norm(r0)>1e-1 && k<=N
 
    ttt = abs(A'*r);        
 
    [~,nn] = max(ttt);
    S(k) = nn;
    col = A(:,S(1:k));
    x1 =(col'*col)\(col'*y);
    r = y - col*x1;
    r0 = r;
    k = k+1;
end
 
x_hat = zeros(N,1);
x_hat(S(1:(k-1))) = x1;
S_hat = S(1:(k-1));
 
end


%%%%%%%%%%%%%%%%%%%% (a) Performance Metrics: %%%%%%%%%%%%%%%%%%%%%%%%%%%
function NMSE = PM(x,x_hat)
NMSE = norm(x-x_hat)/norm(x);
end