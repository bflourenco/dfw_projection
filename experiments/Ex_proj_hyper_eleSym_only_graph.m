clearvars

%% parameters
% problem settings
n=20; deriv_num=10; % Figure 1 in the paper

%n=30; deriv_num=15; % Figure 2 in the paper. Remove the comment marker to
                     % display the plot associated to Figure 2
max_time = 10;
num_points = 10;

rng(1); %fix seed
d = zeros(n,num_points);
poly{1}.type = "symmetric";
poly{1}.deg = n-deriv_num;
poly{1}.index = 1:n;
e = ones(n,1);
for k = 1:num_points
    dk = randn(n,1);
    min_eig_dk = min(eigH(dk,e,poly));
    %Discard test points that are too close to the cone.
    while (min_eig_dk > -1e-4 )
        dk = randn(n,1);
        min_eig_dk = min(eigH(dk,e,poly));
        fprintf("Discarded a test point because it is too close to being feasible\n")
    end
    fprintf("Min eigenvalue of the %dth point: %g \n",k,min_eig_dk)
    d(:,k) = dk;
end

%% measure running time of dual-Frank-Wolfe method using eleSym
fprintf('\n start dual-Frank-Wolfe method using eleSym\n\n')
% set up
clearvars opts
e = ones(n,1);
opts.y0 = zeros(n,1);



for k = 1:num_points
    % measure runtime
    % make problem
    b = -d(:,k);
    step.L = 1;
    step.rule = "Lipschitz";
    opts.step = step;
    grad = @(x) (x-b);
    c = norm(e)*norm(e-d(:,k));
    if dot(e,opts.y0) > c, c = dot(e,opts.y0); end
    
    %stopping criteria
    opts.stop.maxitr = Inf;
    opts.stop.gap_tol = -Inf;
    opts.stop.feas_check = true;
    opts.stop.max_time =max_time;
    [x_FW,y_FW,itr_FW,gap_FW,feas_FW,runtime_FW] = FW_HP_exp(grad,eye(n),zeros(n,1),poly,e,c,opts);
    %calculate objective values
    obj_vals_FW = 0.5*vecnorm(x_FW-d(:,k)).^2;
    FW_obj = obj_vals_FW(feas_FW); 
    

    %% draw graph
    % function value
%     figure(1);
%     %plot(runtime_FW(feas_FW),FW_obj,runtime_AGM,AGM_obj,'--','LineWidth',1);
%     loglog(runtime_FW(feas_FW),FW_obj,"-o","MarkerSize",3); hold on
%     % set legend, xlabel, ylabel, font
%     %legend('FW-EleSym')
%     xlabel('$t(k)$  (unit:second)',"Interpreter","latex")
%     ylabel('$\min_{1\leq i \leq k} f(x_i)$',"Interpreter","latex")
%     %fontsize(gcf,16,"points")
    
    figure(1);
    % FW gap
    %semilogy(runtime_FW(feas_FW),gap_FW(feas_FW),'LineWidth',1);
    %loglog(runtime_FW(feas_FW),gap_FW(feas_FW),'--','LineWidth',1); hold on;
    loglog(runtime_FW(feas_FW),gap_FW(feas_FW)); hold on;
    %yaxisproperties= get(gca, 'YAxis');
    %yaxisproperties.TickLabelInterpreter = 'latex';   
    %yticks([10^-2 10^-1 1 10])
    %yticklabels({'$10^{-2}$','$10^{-1}$','1','10'},"Interpreter","latex");
    grid on;
    %set legend, xlabel, ylabel, font
    xlabel('$t(k)$  (unit:second)',"Interpreter","latex")
    ylabel('FWGap',"Interpreter","latex")
    %fontsize(gcf,16,"points")

    figure(2);
    min_val = min(obj_vals_FW(feas_FW));
    FW_rel = abs(obj_vals_FW(feas_FW)-min_val)/min_val;
    for i = 2:length(runtime_FW(feas_FW))
        FW_obj(i) = min(FW_obj(i-1:i));
        FW_rel(i) = min(FW_rel(i-1:i));
    end
    loglog(runtime_FW(feas_FW),FW_rel); hold on;
    grid on
    %set legend, xlabel, ylabel, font
    xlabel('$t(k)$  (unit:second)',"Interpreter","latex")
    ylabel('$\min_{1\leq i \leq k} (f(x_i) - \hat{f}_{\textrm{opt}})/\hat{f}_{\textrm{opt}}$',"Interpreter","latex")
    %fontsize(gcf,16,"points")
    print(sprintf("objrelplot%d_%d",n,deriv_num),"-deps")

end
%mylinestyles = ["-"; "--"; ":";"-."];

% figure(1);
% %ax = gca;
% %ax.LineStyleOrder = mylinestyles;
% print(sprintf("objplot%d_%d",n,deriv_num),"-deps");
% hold off;

figure(1);
%ax = gca;
%ax.LineStyleOrder = mylinestyles;
print(sprintf("FWgapplot%d_%d",n,deriv_num),"-deps");
hold off;

figure(2);
%ax = gca;
%ax.LineStyleOrder = mylinestyles;
print(sprintf("objrelplot%d_%d",n,deriv_num),"-deps")
hold off;
