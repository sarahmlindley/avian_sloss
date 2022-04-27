%% SLOSS code (3 v 4 patches)
 
% 3 PATCH CODE
clearvars;close all;clc;

% Parameters
mu1    = 0.000000001;    %in-patch mortality for patch 1
mu2    = 0.000000001;    %in-patch mortality for patch 2
mu3    = 0.000000001;    %in-patch mortality for patch 3
phi0   = 0.000000001;    %mortality between patch 1 and 2 
phi1   = 0.000000001;    %mortality between patch 2 and 3
north  = 120;       %days before migration begins
south  = 120;       %days between south and north mig
theta  = 1.25;      %speed power
K      = 30;        %max migration days
a      = 12000;     %linear cost parameter %12000
d      = 0.20;      %migration rate
D      = 500;       %migration distance
b      = 500;       %quadratic cost parameter %500


% Functions 
gamma = @(A) 1./(1+exp(-0.02*A+5))*5; %productivity
C3    = @(A) a*A + b*A^2; %cost per area unit

%5,000,000 to 100,000,000
bud = linspace(10000000,100000000,11);
bb = 1;

% Determining the grid to search over for the manager
mansz = 25;
jjj = 0;
for A = [1,mansz:mansz:500]
    for D1 = 0:mansz:500
        D2 = D-D1;
        jjj = jjj + 1;
        test_vals(jjj,:) = [A,D1,D2];
    end
end

% Grid search over the ecology
% This step generates a set of values to try for tau12, T2, and tau23
    step = K/60; 
    % nn = 0;
    x_gs = [];
    for x1 = 0:step:K
        for x3 = 0:step:K-x1
            x2 = K - x1-x3;
            if(x2>=0)
                x_gs = [x_gs;x1,x2,x3];
            end
        end
    end

% Manager's grid search
for budget = bud
    


% Looping over the manager's grid search
for jj = 1:size(test_vals,1)

% Assigning values to manager's choice variables for each iteration
    A = test_vals(jj,1);
    D1 = test_vals(jj,2);
    D2 = test_vals(jj,3);


% The choice variables for the next three functions are for the ecology
% with the following definitions: x(1) = tau12, x(2) = T2S, x(3) = tau23

% Functions for Survival rate N-S given A,D1,D2
    funNS = @(x) -exp(-mu1*north)*exp(-phi0*exp(theta*D1/x(:,1)))*exp(-mu2*x(:,2))*exp(-phi1*exp(-gamma(A)*x(:,2)+theta*D2/x(:,3)));
    % Log survival rate N-S given A,D1,D2
    logfunNS = @(x) -mu1*north - phi0*exp(theta*D1./x(:,1)) - mu2*x(:,2) - phi1*exp(-gamma(A)*x(:,2) + theta*D2./x(:,3));
    % Log survival rate S-N given A,D1,D2
    logfunSN = @(x) -mu3*south - phi1*exp(theta*D2./x(:,1)) - mu2*x(:,2) - phi0*exp(-gamma(A)*x(:,2) + theta*D1./x(:,3));




    
% Setting constraints for later use of fmincon    
    % Optimization constraints for the ecology
    Aeq = [1,1,1];
    beq = K;
    % Optimization lower and upper bounds for the ecology
    lb = [0,0,0];
    ub = [K,K,K];
    
    
% The north-sourth ecological optimization
    % Determining the log of the survival rate at each of the ecology grid points
    logNS_gs = logfunNS(x_gs);

    % Determining the max of the grid search
    [~,IoptNS] = max(logNS_gs); %tells us where max is in logNS_gs
    xopt_gs_NS(jj,:) = x_gs(IoptNS,:);
    
    % Using the output of the grid search as the initial point for fmincon to
    % find the optimum
    options_fmincon = optimoptions('fmincon','Display','off','OptimalityTolerance',1e-8);
    [xopt_NS(jj,:),fvalNS(jj,1)] = fmincon(@(x) -1*logfunNS(x),xopt_gs_NS(jj,:),[],[],Aeq,beq,lb,ub,[],options_fmincon);
    
% The south-north ecological optimization
    % Determining the log of the survival rate at each of the ecology grid points
    logSN_gs = logfunSN(x_gs);

    % Determining the max of the grid search
    [~,IoptSN] = max(logSN_gs);
    xopt_gs_SN(jj,:) = x_gs(IoptSN,:);
    
    % Using the output of the grid search as the initial point for fmincon to
    % find the optimum
    [xopt_SN(jj,:),fvalSN(jj,1)] = fmincon(@(x) -1*logfunSN(x),xopt_gs_SN(jj,:),[],[],Aeq,beq,lb,ub,[],options_fmincon);

% Determining the total log of the survival rate for the grid search and
% from the outputs of fmincon
    surv_rate(jj,1) = logNS_gs(IoptNS)+logSN_gs(IoptSN);
    surv_rate1(jj,1) = -fvalNS(jj)-fvalSN(jj);

% End manager's grid search loop
end

% generates optimal bird response for every possible value of A, D1, and D2
    % then below manager optimums are calculated given these response rates


% Output
% Manager's optimal values with ecological grid search st to cost
cost = zeros(length(test_vals),1);
    for z = 1:length(test_vals(:,1))
        cost(z,1) = C3(test_vals(z,1));
    end
CSmatrix = [cost,surv_rate,test_vals];
CSmatrix = CSmatrix(CSmatrix(:,1) <= budget , :); %keep all cost =< budget
[~,Imanopt] = max(CSmatrix(:,2)); %find max given budget constraint
manopt(bb,:) = [CSmatrix(Imanopt,1:5)];
% Corresponding ecological values at manager's optimum using grid search
manopt_ecolNS(bb,:) = xopt_gs_NS(Imanopt,:);
manopt_ecolSN(bb,:) = xopt_gs_SN(Imanopt,:);

survivalgs(bb) = exp(-fvalNS(Imanopt))*exp(-fvalSN(Imanopt));

expsurv_rate = exp(surv_rate);
managerall(:,:,bb)=[test_vals,exp(surv_rate),xopt_gs_NS,xopt_gs_SN];

% Manager's optimal values with fmincon
cost1 = zeros(length(test_vals),1);
    for z = 1:length(test_vals(:,1))
        cost1(z,1) = C3(test_vals(z,1));
    end
CSmatrix1 = [cost1,surv_rate1,test_vals];
CSmatrix1 = CSmatrix1(CSmatrix1(:,1) <= budget , :); %keep all cost =< budget
[~,Imanopt1] = max(CSmatrix1(:,2)); %find max given budget constraint
manopt1(bb,:) = [CSmatrix1(Imanopt1,1:5)]; %opt A,D1,D2
% Corresponding ecological values at manager's optimum using fmincon
manopt_ecolNS1(bb,:) = xopt_NS(Imanopt1,:);
manopt_ecolSN1(bb,:) = xopt_SN(Imanopt1,:);

survivalfmc(bb) = exp(-fvalNS(Imanopt1))*exp(-fvalSN(Imanopt1));
bb = bb + 1

end

save("3patch_sloss.mat")
%%
%--------------------------------------------------------------------------
% 4 PATCH CODE

% Parameters
mu1    = 0.000000001;    %in-patch mortality for patch 1
mu2    = 0.000000001;    %in-patch mortality for patch 2
mu3    = 0.000000001;    %in-patch mortality for patch 3
mu4    = 0.000000001;
phi0   = 0.000000001;    %mortality between patch 1 and 2 
phi1   = 0.000000001;    %mortality between patch 2 and 3
phi2   = 0.000000001;
north  = 120;       %days before migration begins
south  = 120;       %days between south and north mig
theta  = 1.25;      %speed power
K      = 30;        %max migration days
a      = 12000;     %linear cost parameter
d      = 0.20;      %migration rate
D      = 500;       %migration distance
b      = 500;       %quadratic cost parameter


% Functions 
gamma  = @(A14) 1/(1+exp(-0.02*A14+5))*5; %productivity
gamma2 = @(A24) 1/(1+exp(-0.02*A24+5))*5; %productivity
C4     = @(A14,A24) a*A14 + b*A14^2 + a*A24 + b*A24^2; %cost per area unit

% Manager's grid search

bud = linspace(1000000,100000000,11);
%bud = [30000000,40000000];
bb4 = 1;

% Determining the grid to search over for the manager
mansz = 25;
jjj = 1;
test_vals4 = zeros();
for A14 = [1,mansz:mansz:250]
     for A24 = [1,mansz:mansz:250]
        for D14 = 0:mansz:500
            for D24 = 0:mansz:500
                D34 = 500 - D24 - D14;
%                 for D34 = 0:mansz:500
%                     if D14 + D24 + D34 == 500
                    if(D34>=0)
                        test_vals4(jjj,1) = A14;
                        test_vals4(jjj,2) = A24;
                        test_vals4(jjj,3) = D14;
                        test_vals4(jjj,4) = D24;
                        test_vals4(jjj,5) = D34;
                        jjj = jjj + 1;
                    end
%                 end
            end
        end
     end
end


load("trial.mat","x_gs4")
% % Grid search over the ecology
% % Generates a set of values to try for tau12, T2, tau23, T3, and tau34 
%     step = K/60; 
%     nn = 1;
%     x_gs4 = zeros();
%     for x14 = 0:step:K
%         for x24 = 0:step:K-x14
%             for x34 = 0:step:K-x14-x24
%                 for x44 = 0:step:K-x14-x24-x34
%                     x54 = K - x14 - x24 - x34 - x44;
%                     if(x54>=0)
%                     x_gs4(nn,1) = x14;
%                     x_gs4(nn,2) = x24;
%                     x_gs4(nn,3) = x34;
%                     x_gs4(nn,4) = x44;
%                     x_gs4(nn,5) = x54;
%                     nn = nn + 1;
%                     if(mod(nn,20000)==0)
%                         nn
%                     end
%                     end
%                 end
%             end
%         end
%     end


for budget = bud


% Looping over the manager's grid search
for jj = 1:size(test_vals4,1)

% Assigning values to manager's choice variables for each iteration
    A14 = test_vals4(jj,1);
    A24 = test_vals4(jj,2);
    D14 = test_vals4(jj,3);
    D24 = test_vals4(jj,4);
    D34 = test_vals4(jj,5);


% The choice variables for the next three functions are for the ecology
% with the following definitions: x(1) = tau12, x(2) = T2S, x(3) = tau23,
% x(4) = T3S, x(5) = tau34, opposite for SN

% Functions for Survival rate N-S given A,D1,D2
    funNS4 = @(x) -exp(-mu1*north)*exp(-phi0*exp(theta*D1/x(:,1)))*exp(-mu2*x(:,2))*exp(-phi1*exp(-gamma(A1)*x(:,2)+theta*D2/x(:,3)))*exp(-mu3*x(:,4))*exp(-phi2*exp(-gamma1(A2)*x(:,4)+theta*D3/x(:,5)));
    % Log survival rate N-S given A,D1,D2
    logfunNS4 = @(x) -mu1*north - phi0*exp(theta*D14./x(:,1)) - mu2*x(:,2) - phi1*exp(-gamma(A14)*x(:,2) + theta*D24./x(:,3)) - mu3*x(:,4) - phi2*exp(-gamma(A24)*x(:,4) + theta*D34./x(:,5));
    % Log survival rate S-N given A,D1,D2
    logfunSN4 = @(x) -mu4*south - phi2*exp(theta*D34./x(:,1)) - mu3*x(:,2) - phi1*exp(-gamma(A24)*x(:,2) + theta*D24./x(:,3)) - mu2*x(:,4) - phi0*exp(-gamma(A14)*x(:,4) + theta*D14./x(:,5));



    
% Setting constraints for later use of fmincon    
    % Optimization constraints for the ecology
    Aeq = [1,1,1,1,1];
    beq = K;
    % Optimization lower and upper bounds for the ecology
    lb = [0,0,0,0,0];
    ub = [K,K,K,K,K];
    
    
% The north-sourth ecological optimization
    % Determining the log of the survival rate at each of the ecology grid points
    logNS_gs4 = logfunNS4(x_gs4);

    % Determining the max of the grid search
    [~,IoptNS4] = max(logNS_gs4); %tells us where max is in logNS_gs
    xopt_gs_NS4(jj,:) = x_gs4(IoptNS4,:);
    
    % Using the output of the grid search as the initial point for fmincon to
    % find the optimum
    options_fmincon = optimoptions('fmincon','Display','off','OptimalityTolerance',1e-8);
    [xopt_NS4(jj,:),fvalNS4(jj,1)] = fmincon(@(x) -1*logfunNS4(x),xopt_gs_NS4(jj,:),[],[],Aeq,beq,lb,ub,[],options_fmincon);

% The south-north ecological optimization
    % Determining the log of the survival rate at each of the ecology grid points
    logSN_gs4 = logfunSN4(x_gs4);

    % Determining the max of the grid search
    [~,IoptSN4] = max(logSN_gs4);
    xopt_gs_SN4(jj,:) = x_gs4(IoptSN4,:);
    
    % Using the output of the grid search as the initial point for fmincon to
    % find the optimum
    [xopt_SN4(jj,:),fvalSN4(jj,1)] = fmincon(@(x) -1*logfunSN4(x),xopt_gs_SN4(jj,:),[],[],Aeq,beq,lb,ub,[],options_fmincon);

% Determining the total log of the survival rate for the grid search and
% from the outputs of fmincon
    surv_rate4(jj,1) = logNS_gs4(IoptNS4)+logSN_gs4(IoptSN4);
    surv_rate14(jj,1) = -fvalNS4(jj)-fvalSN4(jj);

% End manager's grid search loop
end


% generates optimal bird response for every possible value of A, D1, and D2
    % then below manager optimums are calculated given these response rates


% Output
% Manager's optimal values with ecological grid search st to cost
cost4 = zeros(length(test_vals4),1);
    for z = 1:length(test_vals4(:,1))
        cost4(z,1) = C4(test_vals4(z,1),test_vals4(z,2));
    end
CSmatrix4 = [cost4,surv_rate4,test_vals4];
CSmatrix4 = CSmatrix4(CSmatrix4(:,1) <= budget , :); %keep all cost =< budget
[~,Imanopt4] = max(CSmatrix4(:,2)); %find max given budget constraint
manopt4(bb4,:) = [CSmatrix4(Imanopt4,1:7)];
% Corresponding ecological values at manager's optimum using grid search
manopt_ecolNS4(bb4,:) = xopt_gs_NS4(Imanopt4,:);
manopt_ecolSN4(bb4,:) = xopt_gs_SN4(Imanopt4,:);

survival4gs(bb4) = exp(-fvalNS4(Imanopt4)).*exp(-fvalSN4(Imanopt4));

expsurv_rate4 = exp(surv_rate4);
managerall4(:,:,bb4) = [test_vals4,expsurv_rate4,xopt_gs_NS4,xopt_gs_SN4];


% Manager's optimal values with fmincon
cost14 = zeros(length(test_vals4),1);
    for z = 1:length(test_vals4(:,1))
        cost14(z,1) = C4(test_vals4(z,1),test_vals4(z,2));
    end
CSmatrix14 = [cost14,surv_rate14,test_vals4];
CSmatrix14 = CSmatrix14(CSmatrix14(:,1) <= budget , :); %keep all cost =< budget
[~,Imanopt14] = max(CSmatrix14(:,2)); %find max given budget constraint
manopt14 = [CSmatrix14(Imanopt14,:)]; %opt A,D1,D2
% Corresponding ecological values at manager's optimum using fmincon
manopt_ecolNS14(bb4,:) = xopt_NS4(Imanopt14,:);
manopt_ecolSN14(bb4,:) = xopt_SN4(Imanopt14,:);


survival4gs1(bb4) = exp(-fvalNS4(Imanopt14)).*exp(-fvalSN4(Imanopt14));

survival4fmc(bb4) = exp(-fvalNS4(Imanopt14)).*exp(-fvalSN4(Imanopt14));
bb4 = bb4 + 1

end

save("4patch_sloss.mat")

%%
% GRAPHS 

figure;
hold on
plot(bud,exp(manopt(:,2)));
plot(bud,exp(manopt4(:,2)));
xlabel('Budget')
ylabel('Survival Rate')
legend

%%
%% Creating surface plot
load("3patch_sloss.mat")
surf_plot = [managerall(:,1:2),managerall(:,4)];
lngth_side = sqrt(size(surf_plot,1));
lngth_side_1 = length(unique(surf_plot(:,1)));
for i = 1:lngth_side
    X(:,i) = surf_plot((i-1)*lngth_side+1:i*lngth_side,1);
    Y(:,i) = surf_plot((i-1)*lngth_side+1:i*lngth_side,2);
    Z(:,i) = surf_plot((i-1)*lngth_side+1:i*lngth_side,3);
end

% figure;
% contourf(X,Y,Z)
% colorbar;

Z2 = Z;
Z2(Z2<5*surv_rate1(Imanopt1))=5*surv_rate1(Imanopt1);
% figure;
subplot(3,3,gg)
contourf(X,Y,Z2)
colorbar;
title(['$g=$',num2str(gg)],'Interpreter','Latex')
xlabel('$A$','Interpreter','Latex')
ylabel('$D_1$','Interpreter','Latex')
