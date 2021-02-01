spick = [10,15,20];
Mpick = [20,35,50]; % change
Npick = [20,50,100];
f = @PM;
 
%% noiseless
for n = 2:length(Npick)  %% !!!!!!!!
    n=2;
    smax = spick(n);
    Mmax = Mpick(n);
    
    ESRfinal2 = zeros(Mmax,smax);
    ANEfinal2 = zeros(Mmax,smax);
    
    N = Npick(n);
    count = 1;
    Dataall = zeros(smax*Mmax,2);  % make zero matrix !!warning
 
    for M = 1:Mmax
        for s = 1:smax
            ANE = 0;
            ESR = 0;
            for loop = 1:2
                % generate A
                A = randn(M,N);
                B = sqrt(diag(A'*A));
                A = A./B';  % normalize A
                % generate s indices
                S = randi([1 N],s,1);
                S = sort(S); % sort S
                % generate non-zero entries of x
                entries = zeros(s,1);
                for i = 1:s
                    coin = randi([0 1]);
                    if coin == 1
                        entries(i) = 1 + 9*rand;
                    else
                        entries(i) = -1 + (-9)*rand;
                    end
                end
                % generate x
                x = zeros(N,1);
                x(S) = entries;
                % generate y
                y = A(:,S)*entries;
                % noise
                % OMP
                [x_hat,S_hat] = OMP1(A,y);
                % identify S_hat and S, obtain ESR
                if length(S_hat) == s
                    if sum(sort(S_hat)==S)==s
                        ESR = ESR+1;
                    end
                end
           
                % identify ANE(average Normalized Error)
                ANE = ANE + f(x,x_hat);
                
            end
            Dataall(count,:) = [ESR/20,ANE/20]; % 每换一个n就变一次
            count = count+1;
            
            ESRfinal2(M,s)=ESR/20;
            ANEfinal2(M,s)=ANE/20;
        end
        disp(M);    
    end
    
    k = 1;
    ESRfinal = zeros(Mmax,smax);
    ANEfinal = zeros(Mmax,smax);
    for i = 1:Mmax
        ESRfinal(i,:) = Dataall(k:k+smax-1,1);
        ANEfinal(i,:) = Dataall(k:k+smax-1,2);
        k = k+smax;
    end
%     disp(N);
%     figure;
%     imagesc(ESRfinal);
%     colormap gray
%  
%     figure;
%     imagesc(ANEfinal);
%     colormap gray
    
end
%  
%  
%% noise 1
sigmapick = [0.001,0.01]; % sigam^2 = 0.01 and 1
% for n = 1:3  %% !!!!!!!!
%     smax = spick(n);
%     Mmax = Mpick(n);
%     N = Npick(n);
%     for noi = 1:2
%         count = 1;
%         Dataall = zeros(smax*Mmax,1);
%         for M = 1:Mmax
%             noise = sigmapick(noi)*randn(M,1);
%             for s = 1:smax
%                 NEN = 0;
%                 parfor loop = 1:10
%                     % generate A
%                     A = randn(M,N);
%                     B = sqrt(diag(A'*A));
%                     A = A./B';  % normalize A
%                     % generate s indices
%                     S = randi([1 N],s,1);
%                     S = sort(S); % sort S
%                     % generate non-zero entries of x
%                     entries = zeros(s,1);
%                     for i = 1:s
%                         coin = randi([0 1]);
%                         if coin == 1
%                             entries(i) = 1 + 9*rand;
%                         else
%                             entries(i) = -1 + (-9)*rand;
%                         end
%                     end
%                     % generate x
%                     x = zeros(N,1);
%                     x(S) = entries;
%                     % generate y
%                     y = A(:,S)*entries+noise;
%                     % OMP
%                     [x_hat,S_hat] = OMP2(A,y,s);
%                     % judge NEN, normalized error number
%                     if f(x,x_hat)<=1e-3
%                         NEN = NEN + 1;
%                     end
%  
%                 end
%                 Dataall(count) = NEN/10; % 每换一个n就变一次
%                 count = count+1;
%  
%             end
%             disp(M);    
%         end
%         k = 1;
%         NENfinal = zeros(Mmax,smax);
%         for i = 1:Mmax
%             NENfinal(i,:) = Dataall(k:k+smax-1);
%             k = k+smax;
%         end
%         strname = ['noise1N=' num2str(n) 'noi=' num2str(sigmapick(noi)) '.mat'];
%         save(strname,'NENfinal');
%         disp(noi);
%     end
%     
%  
%     disp(N);
%     
%     figure;
%     imagesc(NENfinal);
%     
%     title(strname)
%     colormap gray
%     colorbar
%     
% end
 
%% noise2
 % sigam^2 = 0.01 and 1
for n = 1:3 %% !!!!!!!!
    smax = spick(n);
    Mmax = Mpick(n);
    N = Npick(n);
 
    for noi = 1:2
        Dataall = zeros(smax*Mmax,1);
        count = 1;
        for M = 1:Mmax
            noise = sigmapick(noi)*randn(M,1);
            noiseband = norm(noise);
            for s = 1:smax
                NEN = 0;
                parfor loop = 1:10
                    % generate A
                    A = randn(M,N);
                    A = randn(M,N);
                    B = sqrt(diag(A'*A));
                    % generate s indices
                    S = randi([1 N],s,1);
                    S = sort(S); % sort S
                    % generate non-zero entries of x
                    entries = zeros(s,1);
                    for i = 1:s
                        coin = randi([0 1]);
                        if coin == 1
                            entries(i) = 1 + 9*rand;
                        else
                            entries(i) = -1 + (-9)*rand;
                        end
                    end
                    % generate x
                    x = zeros(N,1);
                    x(S) = entries;
                    % generate y
                    y = A(:,S)*entries+noise;
                    % OMP
                    [x_hat,S_hat] = OMP3(A,y,noiseband);
                    % judge NEN, normalized error number
                    if f(x,x_hat)<=1e-3
                        NEN = NEN + 1;
                    end
 
                end
                Dataall(count) = NEN/10; % 每换一个n就变一次
                count = count+1;
 
            end
            disp(M);    
        end
        k = 1;
        NENfinal = zeros(Mmax,smax);
        for i = 1:Mmax
            NENfinal(i,:) = Dataall(k:k+smax-1);
            k = k+smax;
        end
        strname = ['noise2N=' num2str(n) 'noi=' num2str(sigmapick(noi)) '.mat'];
        save(strname,'NENfinal');
        disp(noi);
        figure;
        imagesc(NENfinal);
        title(strname)
        colormap gray
        colorbar
    end
%     
%     
    disp(N);
 
    
end

% plot
 
    figure;
    imagesc(NENfinal);
    colormap gray
    set(gca,'YDir','normal');
    title({'Noisycase2 Normalized Error; N = 100, $\sigma$ = 0.01'},'Interpreter','latex')
  %  title({'Exact Support Recovery; N = 100'},'Interpreter','tex')
    xlabel({'$sparsity/s$'},'Interpreter','latex');
    ylabel({'$M$'},'Interpreter','latex');
    colorbar
%     axis([0.5 20.5 0.5 50.5]);
%     set(gca,'xtick',[1:2:20.5]);

% Real image case
A = imread('cameraman.tif');
A = A([1:150],[1:150]);
imagesc(A)
Fig = im2double(A);
F_dct = dct2(Fig);
 
%F_dct1
F_dct1 = F_dct;
for i = 1:150
    [~,tindex] = sort(abs(F_dct(:,i)));
    ttt = F_dct(tindex,i);
    [~,Findex] = sort(tindex);
    ttt(1:130) = 0;
    F_dct1(:,i) = ttt(Findex);
end
figure;
imagesc(idct2(F_dct1))
imshow(idct2(F_dct1))
colormap gray
set(gca, 'XTickLabel', [],'XTick',[],'YTickLabel', [],'YTick',[]) 
Fig_n = idct(F_dct1);
sparsity = sum(F_dct1~=0);
 
 
%F_dct2
F_dct2 = F_dct;
F_dct2(abs(F_dct2)<0.01) = 0;
figure;
imagesc(idct(F_dct2))
colormap gray
set(gca, 'XTickLabel', [],'XTick',[],'YTickLabel', [],'YTick',[]) 
 
mratio = [0.2 0.307 0.39 0.5 0.6 0.8];
m = 150*mratio;
for j = 1:length(mratio)
    measure = (1/sqrt(m(j)))*randn(m(j),150);
    x_hat = zeros(150,150);
    for i = 1:150
        yy = measure*F_dct1(:,i);
        [x_hat(:,i),~] = OMP1(measure,yy);
    end
    strratio = ['M/N = ' num2str(mratio(j)) '; sparsity I'];
    subplot(2,3,j);, imagesc(idct(x_hat)), title(strratio)
    set(gca, 'XTickLabel', [],'XTick',[],'YTickLabel', [],'YTick',[]) 
    colormap gray
    
end


function [x_hat, S_hat] = OMP1(A,y)
% S_hat, x_hat are not sorted
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

function [x_hat, S_hat] = OMP2(A,y,s)
% S_hat, x_hat are not sorted
N = size(A,2);
S = zeros(N,1);
r = y;
r0 = 1;
k = 1;
while (k<=s) && ((r0'*r0) > 0.001)
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

function [x_hat, S_hat] = OMP3(A,y,n)
% S_hat, x_hat are not sorted
N = size(A,2);
S = zeros(N,1);
r = y;
r0 = 20;
k = 1;
while r0>n
    ttt = abs(A'*r);  
    [~,nn] = max(ttt);
    S(k) = nn;
    col = A(:,S(1:k));
    x1 = linsolve(col'*col,col'*y);
    r = y - col*x1;
    r0 = norm(r);
    k = k+1;
end
 
x_hat = zeros(N,1);
x_hat(S(1:(k-1))) = x1;
S_hat = S(1:(k-1));
 
end

function y = PM(x,x_hat)
y = norm(x-x_hat)/norm(x);
end

