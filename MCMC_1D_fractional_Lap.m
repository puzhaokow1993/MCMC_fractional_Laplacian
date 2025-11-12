clear all 
close all 
clc 

tic 
% meaningful parameters 
s = 1/2; % order of fractional Laplacian
N = 150; % sample size 
sigma = 0.001; 

% technical parameters for discretization 
M = 50; K=6*M; % number of discretization nodes 
h=6/K; % discretization size 
beta = 0.1; % learning rate
iter = 5000000; % number of MCMC iterations 
J0 = 3; % resolution parameter


% bump function on (-3,-2) in vector form 
phi = zeros(K-1,1); 
for i=1:(M-1)
    phi(i,1) = 10000*exp(1/((-3+i/M+5/2)^2 - (1/4))); 
end

% ---------------------------- FORWARD PROBLEM ----------------------------

% -------------------- discretization of the equation ---------------------
% discretization of fractional Laplacian 
A = zeros(K-1,K-1); 
for i=1:(K-1)
    for k=2:K 
        A(i,i) = A(i,i) + ((k+1)^(1-s)-(k-1)^(1-s))/(k^(1+s)); 
    end
    A(i,i) = A(i,i) + (K^(1-s)-(K-1)^(1-s))/(K^(1+s)) + 2^(1-s) + (1-s)/(s*(K^(2*s)));
    if i<K-1
        A(i,i+1) = -2^(-s); 
    end
    if i>1 
        A(i,i-1) = -2^(-s); 
    end
end
for i=1:(K-1)
    for j=1:(K-1)
        if abs(i-j)>1 
            A(i,j) = -((abs(j-i)+1)^(1-s)-(abs(j-i)-1)^(1-s))/(2*abs(j-i)^(1+s)); 
        end
    end
end
A = A*((4^s)*s*gamma((1+2*s)/2))/((sqrt(pi)*gamma(1-s))*(1-s)*(h^(2*s))); 

fground = eye(K-1); % the diagonal elements are discretization of the function f (initial guess) 
f = eye(K-1); % the diagonal elements are discretization of the function f (to be reconstructed)  
for i=1:(K-1)
    if (i>(5*M)/2) && (i<(7*M)/2) 
        f(i,i) = 1+100*exp(1/((-3+i*h)^2 - (1/4))); 
    end 
end

Bground = A(2*M+1:4*M-1,2*M+1:4*M-1) + fground(2*M+1:4*M-1,2*M+1:4*M-1); 
F1 = -A*phi; 
F = F1(2*M+1:4*M-1,1); 
v0ground = Bground\F; % solve equation
vground = zeros(K-1,1); 
vground(2*M+1:4*M-1) = v0ground; 
B = A(2*M+1:4*M-1,2*M+1:4*M-1) + f(2*M+1:4*M-1,2*M+1:4*M-1); 
% F1 = -A*phi; 
% F = F1(2*M+1:4*M-1,1); 
v0 = B\F; % solve fractional equation
v = zeros(K-1,1); 
v(2*M+1:4*M-1) = v0; 

axis = zeros(K-1,1); 
for i=1:(K-1)
    axis(i,1) = -3 + i*h; 
end


figure
plot(axis,vground,'r-','LineWidth',2) 
hold on 
plot(axis,v,'b-','LineWidth',2) 
legend({'$f(x)=1$','$f(x)=1+100\exp(\frac{1}{x^2 - \frac{1}{4}})\chi_{(-\frac{1}{2},\frac{1}{2})}(x)$'},'Interpreter','latex','FontSize',16,'Location','northoutside')
box off 
xlim([-1,1]); 

Gfullground = A*(vground+phi); 
forwardground = [Gfullground(1:2*M-1,1),Gfullground(4*M+1:K-1,1)]; 

Gfull = A*(v+phi); 
forward = [Gfull(1:2*M-1,1);Gfull(4*M+1:K-1,1)]; % exterior measurement which used to reconstruct f

% figure 
% plot(axis(4*M+1:K-1,1),forwardground,'r-','LineWidth',2) 
% hold on 
% plot(axis(4*M+1:K-1,1),forward,'b-','LineWidth',2) 
% legend({'$f(x)=1$','$f(x)=1+100\exp(\frac{1}{x^2 - \frac{1}{4}})\chi_{(-\frac{1}{2},\frac{1}{2})}(x)$'},'Interpreter','latex','FontSize',16,'Location','northoutside')
% box off 

% ---------------------------- INVERSE PROBLEM ----------------------------
% Goal: reconstruct the potential 'f' from the exterior measurement 'forward' 

% simulating the sampling experiment 
randompoint = ceil((4*M-2)*rand(N,1)); 
sample = forward(randompoint); % sampling experiment 

% initial guess: constant function f(x)=1 
sampleinit = forwardground(randompoint); % using the forward operator as mentioned above, and we can compute the measurement on the chosen random points 

array_log_likelihood = zeros(1,iter);
f_seq = cell(1,iter+1); % create an empty cell 
f_seq{1} = fground; % loading initial guess  
keys = cell(1,iter+1); % create an empty cell 
keys{1} = 'accept'; 

for tau=1:iter
    % loading previous iteration 
    f_initial = f_seq{tau}; 
    
    % propose f 
    f_seq_next = f_initial; 
    for r=-2^(J0):1:(2^(J0)-1)
        temp = randn; 
        for i=1:(K-1) 
            f_seq_next(i,i) = 1 + sqrt(1 - beta^2)*(f_seq_next(i,i)-1) + beta.*temp.*(2*(2^J0)*(-3+i*h)-r>0).*(2*(2^J0)*(-3+i*h)-r<=1); 
        end 
    end

    % compute log-likelihood if previous f accept 
    if strcmp(keys{tau}, 'accept')
        B = A(2*M+1:4*M-1,2*M+1:4*M-1) + f_initial(2*M+1:4*M-1,2*M+1:4*M-1); 
        % F1 = -A*phi; 
        % F = F1(2*M+1:4*M-1,1); 
        v0 = B\F; % solve equation in high precision 
        v = zeros(K-1,1); 
        v(2*M+1:4*M-1) = v0; 
        Gfull = A*(v+phi); 
        forward = [Gfull(1:2*M-1,1),Gfull(4*M+1:K-1,1)]; % exterior measurement which used to reconstruct f
        forward_sample = forward(randompoint); 
        current_log_likelihood = norm(sample-forward_sample)/(-sqrt(2)*sigma); 
    end
    array_log_likelihood(1,tau) = current_log_likelihood; 
    
    % compute log-likelihood for proposed f 
    B = A(2*M+1:4*M-1,2*M+1:4*M-1) + f_seq_next(2*M+1:4*M-1,2*M+1:4*M-1); 
    % F1 = -A*phi; 
    % F = F1(2*M+1:4*M-1,1); 
    v0 = B\F; % solve equation 

    v = zeros(K-1,1); 
    v(2*M+1:4*M-1) = v0; 
    Gfull = A*(v+phi); 
    forward = [Gfull(1:2*M-1,1),Gfull(4*M+1:K-1,1)]; % exterior measurement which used to reconstruct f
    forward_sample = forward(randompoint); 
    proposed_log_likelihood = (norm(sample-forward_sample))/(-sqrt(2)*sigma); 
    alpha = proposed_log_likelihood - current_log_likelihood; 

    if (alpha > 0)  
        f_seq{tau+1} = f_seq_next; % accept the proposal f_seq_next if there is an improvement
        keys{tau+1} = 'accept'; 
    else
        f_seq{tau+1} = f_seq{tau}; % reject the proposal f_seq_next  
        keys{tau+1} = 'reject'; 
    end
    disp(['Iteration: ', num2str(tau),'/', num2str(iter)]) % show the number of current iteration 
    disp(['current log likelihood: ', num2str(current_log_likelihood)]) % show the number of current iteration 
    % if mod(tau,50000) == 0 
    %     save('backup') 
    % end
    f_seq{tau} = diag(f_seq{tau}); % saving only the diagonal elements, since all nondiagonal elements are zero 
end 
f_seq{iter+1} = diag(f_seq{iter+1}); % saving only the diagonal elements, since all nondiagonal elements are zero 

key = 'accept'; % key to extract 
idx = strcmp(keys,key); % find indices of matching keys 
f_accepted = f_seq(idx); % extract corresponding values 

% Define the burn-in period (number of samples to discard) 
burnIn = floor(length(f_accepted)/2); 
f_accepted_remain = f_accepted(burnIn+1:end); % extract the values after the burn-in 

sum_f = f_accepted_remain{1}; 
for h=2:length(f_accepted_remain) 
    sum_f = sum_f + f_accepted_remain{h}; 
end
burn_in_mean_f = sum_f/length(f_accepted_remain); 
elapsedtime = toc; 

% save all workspace variables in the file experiment1.mat 
save('experiment1', '-v7.3') 


% Plot of result 
figure 
plot(axis,burn_in_mean_f,'r-','LineWidth',2) 
hold on 
plot(axis,diag(f),'b-','LineWidth',2) 
legend({'burn-in mean $f_{\rm burn}$','true potential $f$'},'Interpreter','latex','FontSize',16,'Location','northoutside')
box off 
xlim([-1,1]); 

% Plot log-likelihood of each iterations (log scale on x axis)  
figure; 
iter_array = 1:iter; 
semilogy(iter_array,array_log_likelihood,'-o')
xlim([10001,iter])