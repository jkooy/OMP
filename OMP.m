%%%%%%%%%%%%%%%%%%%%%%%%%  ECE269 project  XUEHAI HE %%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;
times=2000;

%%%%%%%%%%%%%%%%%%%%%%%%% (c) Noiseless case: (n = 0)%%%%%%%%%%%%%%%%%%%%%%%%%%


[ESR_matrix, aNE_matrix]= Monte_Carlo_runs(0, times, 1);


%%%%%%%%%%%%%%%%%%%%%%%%% (d) Noisy case: (n != 0) %%%%%%%%%%%%%%%%%%%%%%%%

    sigma_list = [1e-3,1e-2]; 
    % two values of s (one small and one large)


%% a sparsity s is known
    
    for sigma= 1:2    
        [ESR_matrix, aNE_matrix]= Monte_Carlo_runs(sigma_list(sigma), times, 1);
    end

    
%% b sparsity s is NOT known
    for sigma= 1:2    
        [ESR_matrix, aNE_matrix]= Monte_Carlo_runs(sigma_list(sigma), times, 0);
    end

    
    
%% c test OMP on real images
lena_3 = imread('lena.jpg');
%to single channel
lena = lena_3(:,:,1);

%look at the single channel image
% imagesc(lena)

%normalization to [0,1] with double
x = im2double(lena);

%We find the DCT matrix for which the image x has a sparse representation
DCT=dctmtx(size(lena_3,1));

%create the sparse input image
s=DCT\x;
s(abs(s)<0.2) = 0;

figure;
%let's see the sparse input image
x_sparse=DCT*s;
imagesc(x_sparse);
colormap gray;
sparsity = sum(sum(s~=0))/numel(s);
 

figure; 
measurement_rate = [0.2 0.31 0.39 0.5 0.6 0.65 0.7 0.75 0.8];
m = 200*measurement_rate;
for j = 1:length(measurement_rate)
    %corrupt the image with a random matrix A
    A = randn(m(j),200);
    s_hat = zeros(200,200);
    B=A*DCT;
    for i = 1:200
        y_o = A*x_sparse(:,i);
        
        %do not add noise
         
        [s_hat(:,i),~] = OMP1(B,y_o);
    end
    subplot(3,3,j);
    imagesc(DCT*s_hat);
%     imagesc(idct2(x_hat));
    title(['MR= ' num2str(measurement_rate(j)) ', Sparsity=',num2str(sparsity)]);
    colormap gray
    
end




function [ESR_matrix, aNE_matrix]= Monte_Carlo_runs(sigma, times, sparsity_known)
    for n = 1:4
    %% d
        % Generate different phase transition plots for the following values of N: 20, 50 and 100.
        N_list=[20,50,80,100];
        % We generate for N: 20, 50, 80 and 100.
        %the average Normalized Error should be repeating step 1) to step 3) 2000 times and averaging the results over these 2000 Monte Carlo runs
        measurement_rate = [1,0.8,0.6,0.5];
        M_list=N_list.*measurement_rate;
        Sparsity_max = [0.5,0.24,0.2,0.18];
        S_list = N_list.*Sparsity_max;
        
        
        N=N_list(n);
        M_max=M_list(n);
        S_max=S_list(n);
%         M_max=35;
        % Smax is chosen to be a reasonably large integer which is smaller than N
%         S_max=0.3*N;
%         S_max=15;


        ESR_matrix = zeros(M_max,S_max);
        aNE_matrix = zeros(M_max,S_max);
        NE_matrix = zeros(M_max,S_max);

        for M = 1:M_max
            noise = randn(M,1);
            normnoise = norm(sigma*noise);
            
            
            for s = 1:S_max 

                ESR=0;
                aNE=0;
                NE=0;
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
                    for index = 1:s
                        U = randi([0 1]);
                        if U==1
                            e(index) = 1+9*rand;
%                             in the range [􀀀10,􀀀1] [ [1, 10].
                        elseif U==0
                            e(index) = -1-9*rand;
                        end
                    end
                    % generate the original signal x
                    x = zeros(N,1);
                    x(S) = e;
                    %sorted signal                   
                    SNR=norm(x)/sigma;
                    %calculate signal to noise ratio
                    SNR_dB=10*log(SNR);
    %% c 
                    
                    % Noiseless case(n = 0)
                    
                    y = A(:,S)*e+sigma*noise;
                    %Noisy case: (n != 0)
                    
                    
                    if sigma==0 && sparsity_known==1
                        [x_hat,S_hat] = OMP1(A,y);

                   %Compute normalized error and exact support recovery for each iteration.
                        if length(S_hat) == s
                            %judge whether they are the same
                            if sum(sort(S_hat)==S)==s
                                ESR = ESR+1;
                            end
                        end
                        % identify ANE(average Normalized Error)
                        aNE = aNE + PM(x,x_hat);   
                          
                    elseif sigma~=0 && sparsity_known==1
                        x_hat = OMP_s(A,y,s);
                        % success is defined as the event that the Normalized Error is less than 10e-3
                        if PM(x,x_hat)<=1e-3
                            NE = NE + 1;
                        %Success add 1.
                        end
                        
                    elseif sigma~=0 && sparsity_known==0
                        x_hat = OMP_n2(A,y,normnoise);
                        %success is defined as the event that the Normalized Error is less than 10e-3.
                        if PM(x,x_hat)<=1e-3
                            NE = NE + 1;
                        %Success add 1.
                        end
                        
                    end 
                        
                end
                ESR_matrix(M,s)=ESR/times;
                aNE_matrix(M,s)=aNE/times;
                NE_matrix(M,s)=NE/times;
            end
        end
        
        if sigma==0
            if n==1
               figure;
            end
            %noiseless phase transition
            subplot(2,2,n);
            imagesc(ESR_matrix);
            colormap gray;
            title(['noiseless phase transition (ESR)']);
            xlabel(['Smax (N=', num2str(N), ', sigma=' ,num2str(sigma),')']);
            ylabel('M');
            colorbar;
            set(gca,'YDir','normal');

%             if n==1
%                figure;
%             end
%        %Regenerate phase transition plots for average Normalized Error(instead of probability of successful recovery）
%             subplot(2,2,n);
%             imagesc(aNE_matrix);
%             colormap gray;
%             title('noiseless phase transition (AVE)');
%             xlabel(['Smax (N=', num2str(N), ', sigma=' ,num2str(sigma),')']);
%             ylabel('M');
%             colorbar;
%             set(gca,'YDir','normal');
            if n==1
                saveas(gcf,['ESR_noiseless_Smax (N=', num2str(N), '_sigma=' ,num2str(sigma),').png']);
            end
        else
%             figure;
            if n==1
               figure;
            end
            subplot(2,2,n);
            imagesc(NE_matrix);
            colormap gray;
            title('noisy phase transition (Normalized Error)');
            xlabel(['Smax (N=', num2str(N), ', sigma=' ,num2str(sigma),')']);
            ylabel('M');
            colorbar;
            set(gca,'YDir','normal');
            if n==1
                if sparsity_known==1
                    saveas(gcf,['S_noisy_Smax (N=', num2str(N), '_sigma=' ,num2str(sigma),').png']);
                elseif sparsity_known==0
                    saveas(gcf,['n_noisy_Smax (N=', num2str(N), '_sigma=' ,num2str(sigma),').png']);
                end  
            end
        end
    end
end



function [x_hat, S_hat] = OMP1(A,y)
N = size(A,2);
Support = zeros(N,1);
inde = 1;
%initialize the residue with a large number
residue_0 = 10000;

while (norm(residue_0)>1e-1) && (inde<=N)
    product = abs(A'*residue);        
 
    [~,n] = max(product);
    Support(inde) = n;
    A_col = A(:,Support(1:inde));
    %update A column
    x1 =(A_col'*A_col)\(A_col'*y);
    residue_0= y - A_col*x1;
%     calculate the residue
    inde = inde+1;
end
 
x_hat = zeros(N,1);
x_hat(Support(1:(inde-1))) = x1;
S_hat = Support(1:(inde-1));
 
end


% function x_hat = OMP_s(A,y,s)
% % Assume that sparsity s is known
% %initialize support, index, residual 
% residue = y;
% N = size(A,2);
% Support = zeros(N,1);
% residue_0 = 1;
% inde = 1;
% while (inde<=s) && ((residue_0'*residue_0) > 0.001)
%     %terminate the algorithm after first s columns of A are selected
%     product = abs(A'*residue);
%     [~,n] = max(product);
%     Support(inde) = n;
%     A_col = A(:,Support(1:inde));
%     %update A column
%     x =(A_col'*A_col)\(A_col'*y);
%     residue = y - A_col*x;
% %     calculate the residue
%     residue_0 = residue;
%     inde = inde+1;
% end
% x_hat = zeros(N,1);
% x_hat(Support(1:(inde-1))) = x;
%  
% end
% % 
% % 
% function x_hat = OMP_n2(A,y,norm)
% % the sparsity s is NOT known, but ||n||2 is known
% %initialize support, index, residual
% residue = y;
% N = size(A,2);
% Support = zeros(N,1);
% %initialize the residue with a large number
% residue_0 = 10000;
% index = 1;
% while residue_0>norm
% %   judge the norm of residue, if smaller than the norm of the noise, then break
%     %stop the OMP iterations once ||y-Ax(k)||^2<=||n||2.
%     product = abs(A'*residue);  
%     [~,n] = max(product);
%     Support(index) = n;
%     A_col = A(:,Support(1:index));
%     %update A column
%     x= linsolve(A_col'*A_col,A_col'*y);
%     residue = y - A_col*x;
% %     calculate the residue
%     residue_0 = norm(residue);
% %       judge the norm of residue
%     index = index+1;
% end
% x_hat = zeros(N,1);
% x_hat(Support(1:(index-1))) = x;
%  
% end

%%%%%%%%%%%%%%%%%%%% (a) Performance Metrics: %%%%%%%%%%%%%%%%%%%%%%%%%%%
function NMSE = PM(x,x_hat)
NMSE = norm(x-x_hat)/norm(x);
end