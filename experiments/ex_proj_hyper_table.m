function[] = ex_proj_hyper_table(n,deriv_num,num_points,DDS_tol,ratio_list,save_to_file)
%This is the script for the experiments of projection onto hyperbolicity
%cones corresponding to the derivative relaxations of the nonnegative
%orthant.
%
% If save_to_file is true, a mat file is created with all the data 
% generated during the experiments. Furthermore, a csv file is created 
% with the experiment statistics. 
%
% Example:ex_proj_hyper_table(30,27,10,1e-8,[1e-1,5e-2,1e-2,5e-3,1e-3,5e-4,1e-4,5e-5,1e-5],false)

% Generate points to project
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% measure running time of DDS
fprintf('Start DDS experiments \n\n')

% set up
runtime_DDS = zeros(num_points,1); 
x_DDS = cell(1,num_points);
obj_vals_DDS = cell(1,num_points);
itr_DDS = cell(1,num_points);

% measure running time
for k = 1:num_points
    fprintf('%d-th point\n\n',k)
    %%%%% build the problem %%%%%
    u = d(:,k);

    % add a SOCP constraint
    cons{1,1}='SOCP';
    cons{1,2}=[n];
    A{1,1}=[1,zeros(1,n);
            zeros(n,1),eye(n)];
    b{1,1}=[0;-u];
  
    %make poly for hyperbolicity cone constraint in DDS
    combi = nchoosek(1:n,n-deriv_num);
    row = size(combi,1);
    poly = zeros(row,n+1);
    for i = 1:row
        for j = combi(i,:)
            poly(i,j) = 1;
        end
    end
    poly(:,n+1) = ones(row,1);
    % DDS paper suggests using sparse matrices (section 11.1 in their MPC paper), 
    % so we sparsify the matrix if it has 25% or less nonnzeros
    %
    if nnz(poly)/numel(poly) <= 0.25
        poly = sparse(poly);
    end
    % adding HB constraint 
    cons{2,1}='HB';
    cons{2,2}=[n];
    cons{2,3}=poly;
    cons{2,4}='monomial';
    cons{2,5}=[ones(n,1)];
    A{2,1}=[zeros(n,1),eye(n)];
    b{2,1}=[zeros(n,1)];
 
    % set cost vector c
    c=[1;zeros(n,1)];
    
    OPTIONS.tol = DDS_tol;
    
    [x,~,info]=DDS(c,A,b,cons,OPTIONS);
    x_DDS{k} = x(2:n+1);
    runtime_DDS(k) = info.time;
    itr_DDS{k}(end+1) = info.iter;
    obj_vals_DDS{k} = 0.5*vecnorm(x_DDS{k}-u).^2;
end

%% measure running time of dual-Frank-Wolfe method
fprintf('\n start dual-Frank-Wolfe method\n\n')
% set up
clearvars poly 
e = ones(n,1);

% make poly for hyperbolicity cone constraint
combi = nchoosek(1:n,n-deriv_num);
row = size(combi,1);
poly{1}.type = "matrix";
poly{1}.mat = zeros(row,n+1);
for i = 1:row
    for j = combi(i,:)
        poly{1}.mat(i,j) = 1;
    end
end
poly{1}.mat(:,n+1) = ones(row,1);

%%measure runtime
%prepare data strorages for all datasets d
x_FW = cell(1,num_points);
y_FW = cell(1,num_points);
itr_FW = zeros(1,num_points);
gap_FW = cell(1,num_points);
feas_FW = cell(1,num_points);
runtime_FW = cell(1,num_points);
obj_vals_FW = cell(1,num_points);

%measure runtime for each d
for k = 1:num_points
    fprintf('\n %d-th point\n\n',k)
    clear("opts");
    opts.y0 = zeros(n,1);
    %make problem
    b = -d(:,k);
    step.rule = "Lipschitz";
    step.L = 1;
    opts.step = step;
    grad = @(x) (x-b);
    %c = norm(e)*sqrt(norm(e-d(:,k)).^2+norm(e).^2);
    c = norm(e)*norm(e-d(:,k));
    if dot(e,opts.y0) > c, c = dot(e,opts.y0); end
    
    %opts.zero_eps = 1e-1;
    %opts.verbose_freq = 10;

    %stopping criteria
    opts.stop.maxitr = Inf;
    opts.stop.feas_check = true;

    %It is necessary to stop the algorithm when the FW gap is too small or is
    %zero, since this may lead to numerical problems.
    %Still, we want to algorithm to keep going until we reach the desired 
    %precision.
    %Setting it to 1e-12 ensures that it will rarely (or never) stop 
    %because of the FW gap stopping criterion.
    opts.stop.gap_tol = 1e-12;

    %in order to save time, we stop the experiment if an iterate is 
    %already within min(ratio_list) of the objective value obtained by DDS.
    opts.stop.f = @(x) (0.5*norm(x-d(:,k)).^2);
    opts.stop.obj_thresh = min(obj_vals_DDS{k})*(1+min(ratio_list));
    opts.stop.max_time = runtime_DDS(k); 
  
    [x_FW{k},y_FW{k},itr_FW(k),gap_FW{k},feas_FW{k},runtime_FW{k}] = FW_HP_exp(grad,eye(n),zeros(n,1),poly,e,c,opts);
    
    %calculate objective values
    obj_vals_FW{k} = 0.5*vecnorm(x_FW{k}-d(:,k)).^2;
end

%% measure running time of dual-Frank-Wolfe method using eleSym
fprintf('\n start dual-Frank-Wolfe method using eleSym\n\n')
% remake poly using "symmetric" option
clearvars poly
poly{1}.type = "symmetric";
poly{1}.deg = n-deriv_num;
poly{1}.index = 1:n;

%%measure runtime
%prepare data strorages for all datasets d
x_FW_ele = cell(1,num_points);
y_FW_ele = cell(1,num_points);
itr_FW_ele = zeros(1,num_points);
gap_FW_ele = cell(1,num_points);
feas_FW_ele = cell(1,num_points);
runtime_FW_ele = cell(1,num_points);
obj_vals_FW_ele = cell(1,num_points);

%measure runtime for each d
for k = 1:num_points
    fprintf('\n %d-th point\n\n',k)
    clear("opts");
    opts.y0 = zeros(n,1);
    %make problem
    b = -d(:,k);
    step.L = 1;
    step.rule = "Lipschitz";
    opts.step = step;
    grad = @(x) (x-b);
    c = norm(e)*norm(e-d(:,k));
    if dot(e,opts.y0) > c, c = dot(e,opts.y0); end
    
    %opts.zero_eps = 1e-1;
    %opts.verbose_freq = 10;

    %stopping criteria
    opts.stop.maxitr = Inf;
    
    %It is necessary to stop the algorithm when the FW gap is too small or is
    %zero, since this may lead to numerical problems.
    %Still, we want to algorithm to keep going until we reach the desired 
    %precision.
    %Setting it to 1e-12 ensures that it will rarely (or never) stop 
    %because of the FW gap stopping criterion.
    opts.stop.gap_tol = 1e-12;
   
    opts.stop.feas_check = true;
    %in order to save time, we stop the experiment if an iterate is 
    %already within min(ratio_list) of the objective value obtained by DDS.
    opts.stop.f = @(x) (0.5*norm(x-d(:,k)).^2);
    opts.stop.obj_thresh = min(obj_vals_DDS{k})*(1+min(ratio_list));
    opts.stop.max_time = runtime_DDS(k);
  
    [x_FW_ele{k},y_FW_ele{k},itr_FW_ele(k),gap_FW_ele{k},feas_FW_ele{k},runtime_FW_ele{k}] = FW_HP_exp(grad,eye(n),zeros(n,1),poly,e,c,opts);
    
    %calculate objective values
    obj_vals_FW_ele{k} = 0.5*vecnorm(x_FW_ele{k}-d(:,k)).^2;
end


%% parameters for Renegar's method %%
L = 1;

%% measure running time of Accelerated gradient method
fprintf("\n start Renegar's method\n\n")
%setup
clearvars poly A c e opts

%make poly for hyperbolicity cone constraint in AGM
combi = nchoosek(1:n,n-deriv_num);
row = size(combi,1);
poly{1}.type = "matrix";
poly{1}.mat = zeros(row,n);
for i = 1:row
    for j = combi(i,:)
        poly{1}.mat(i,j) = 1;
    end
end
poly{1}.mat = [zeros(row,n+1),poly{1}.mat,ones(row,1)];

poly{2}.type = "oracle";
poly{2}.deg = 2;
poly{2}.p = @(x) (x(1).^2-sum(x(2:n+1).^2));
poly{2}.grad = @(x) ([2*x(1);-2*x(2:n+1);zeros(n,1)]);

%measure runtime
x_AGM = cell(1,num_points);
itr_AGM = zeros(1,num_points);
num_grad_cal = zeros(1,num_points);
runtime_AGM = cell(1,num_points);
obj_vals_AGM = cell(1,num_points);
for k = 1:num_points
    fprintf('\n %d-th point\n\n',k)
    clear("opts");
    % make problem
    c = [1;zeros(2*n,1)];
    A = [zeros(n,1),-eye(n),eye(n)];
    e = [norm(ones(n,1)-d(:,k))+1;ones(n,1)-d(:,k);ones(n,1)];
    v0 = [norm(ones(n,1)-d(:,k));ones(n,1)-d(:,k);ones(n,1)];

    % stopping criteria
    opts.max_itr = Inf;
    opts.max_num = Inf;
    opts.f_threshold.f = @(x) (0.5*norm(x(2:n+1)-d(:,k)).^2);
    opts.f_threshold.threshold = min(obj_vals_DDS{k})*(1+min(ratio_list));
    opts.max_time = runtime_DDS(k);%%%%%

    % measure running time
    [x_AGM{k},~,~,itr_AGM(k),num_grad_cal(k),runtime_AGM{k}] = AGM_HP(c,A,v0,poly,e,L,opts);
    x_AGM{k} = x_AGM{k}(n+2:end,:);
    obj_vals_AGM{k} = 0.5*vecnorm(x_AGM{k}-d(:,k)).^2;
end

%% measure running time of Accelerated gradient method using eleSym
fprintf("\n start Renegar's method using eleSym\n\n")
%setup
clearvars poly

%make poly for hyperbolicity cone constraint in AGM
poly{1}.type = "symmetric";
poly{1}.deg = n-deriv_num;
poly{1}.index = [n+2:2*n+1];

poly{2}.type = "oracle";
poly{2}.deg = 2;
poly{2}.p = @(x) (x(1).^2-sum(x(2:n+1).^2));
poly{2}.grad = @(x) ([2*x(1);-2*x(2:n+1);zeros(n,1)]);

%measure runtime
x_AGM_ele = cell(1,num_points);
itr_AGM_ele = zeros(1,num_points);
num_grad_cal_ele = zeros(1,num_points);
runtime_AGM_ele = cell(1,num_points);
obj_vals_AGM_ele = cell(1,num_points);
for k = 1:num_points
    fprintf('\n %d-th point\n\n',k)
    clear("opts");
    % build problem data
    c = [1;zeros(2*n,1)];
    A = [zeros(n,1),-eye(n),eye(n)];
    e = [norm(ones(n,1)-d(:,k))+1;ones(n,1)-d(:,k);ones(n,1)];
    v0 = [norm(ones(n,1)-d(:,k));ones(n,1)-d(:,k);ones(n,1)];

    % set up stopping criteria
    opts.max_itr = Inf;
    opts.max_num = Inf;
    opts.f_threshold.f = @(x) (0.5*norm(x(2:n+1)-d(:,k)).^2);
    opts.f_threshold.threshold = min(obj_vals_DDS{k})*(1+min(ratio_list));
    opts.max_time = runtime_DDS(k);%%%%%

    % measure running time
    [x_AGM_ele{k},~,~,itr_AGM_ele(k),num_grad_cal_ele(k),runtime_AGM_ele{k}] = AGM_HP(c,A,v0,poly,e,L,opts);
    x_AGM_ele{k} = x_AGM_ele{k}(n+2:end,:);
    obj_vals_AGM_ele{k} = 0.5*vecnorm(x_AGM_ele{k}-d(:,k)).^2;
end

%% calculate mean time to achieve some relative errors
t_FW = NaN(length(ratio_list),num_points);
t_FW_ele = NaN(length(ratio_list),num_points);
t_AGM = NaN(length(ratio_list),num_points);
t_AGM_ele = NaN(length(ratio_list),num_points);
success_FW = zeros(length(ratio_list),1);
success_FW_ele = zeros(length(ratio_list),1);
success_AGM = zeros(length(ratio_list),1);
success_AGM_ele = zeros(length(ratio_list),1);
for i = 1:length(ratio_list)
    for k = 1:num_points
        if ~isempty(runtime_FW{k}( (obj_vals_FW{k}-obj_vals_DDS{k})<=obj_vals_DDS{k}*ratio_list(i) & feas_FW{k} ))
            t_FW(i,k) = min(runtime_FW{k}((obj_vals_FW{k}-obj_vals_DDS{k})<=obj_vals_DDS{k}*ratio_list(i) & feas_FW{k}) );
            success_FW(i) = success_FW(i)+1;
        end
        if ~isempty(runtime_FW_ele{k}((obj_vals_FW_ele{k}-obj_vals_DDS{k})<=obj_vals_DDS{k}*ratio_list(i) & feas_FW_ele{k}))
            t_FW_ele(i,k) = min(runtime_FW_ele{k}((obj_vals_FW_ele{k}-obj_vals_DDS{k})<=obj_vals_DDS{k}*ratio_list(i) & feas_FW_ele{k} ));
            success_FW_ele(i) = success_FW_ele(i)+1;
        end
        if ~isempty(runtime_AGM{k}((obj_vals_AGM{k}-obj_vals_DDS{k})<=obj_vals_DDS{k}*ratio_list(i)))
            t_AGM(i,k) = min(runtime_AGM{k}((obj_vals_AGM{k}-obj_vals_DDS{k})<=obj_vals_DDS{k}*ratio_list(i)));
            success_AGM(i) = success_AGM(i)+1;
        end
        if ~isempty(runtime_AGM_ele{k}((obj_vals_AGM_ele{k}-obj_vals_DDS{k})<=obj_vals_DDS{k}*ratio_list(i)))
            t_AGM_ele(i,k) = min(runtime_AGM_ele{k}((obj_vals_AGM_ele{k}-obj_vals_DDS{k})<=obj_vals_DDS{k}*ratio_list(i)));
            success_AGM_ele(i) = success_AGM_ele(i)+1;
        end
    end
end

% display results
fprintf("\n Average time to achieve relative errors with respect to DDS\n\n")
fprintf("relative error;\tFW;std;success;\t\tFW_ele;std;success;\tAGM;std;success;\tAGM_ele;std;success\n")
for i = 1:length(ratio_list)
    fprintf("%.2e;\t%4.2f;%4.2f;%4.1f;\t%4.2f;%4.2f;%4.1f;\t%4.2f;%4.2f;%4.1f;\t%4.2f;%4.2f;%4.1f\n",ratio_list(i),...
        mean(t_FW(i,:)'./runtime_DDS*100,"omitnan"),std(t_FW(i,:)'./runtime_DDS*100,"omitnan"),success_FW(i)/num_points*100,...
        mean(t_FW_ele(i,:)'./runtime_DDS*100,"omitnan"),std(t_FW_ele(i,:)'./runtime_DDS*100,"omitnan"),success_FW_ele(i)/num_points*100,...
        mean(t_AGM(i,:)'./runtime_DDS*100,"omitnan"),std(t_AGM(i,:)'./runtime_DDS*100,"omitnan"),success_AGM(i)/num_points*100,...
        mean(t_AGM_ele(i,:)'./runtime_DDS*100,"omitnan"),std(t_AGM_ele(i,:)'./runtime_DDS*100,"omitnan"),success_AGM_ele(i)/num_points*100)
end

if save_to_file
    filename = "proj_hyper_"+sprintf("n%d_",n)+sprintf("d%d",deriv_num)+sprintf("_%d",num_points);
    if (DDS_tol < 1e-3)
        filename = filename + sprintf("_tol_high");
    else
        filename = filename + sprintf("_tol_low");
    end
    datetime.setDefaultFormats('default','yyyy_MM_dd_hh_mm_ss');
    filename = filename + "_" + string(datetime("now"));
    save(filename + ".mat");
    datetime.setDefaultFormats('reset');
    
    fileID = fopen(filename +"stats" +".csv",'w');
    %fprintf(fileID,"n;d;tol;points\n");
    %fprintf(fileID,"%d;%d;%.1e;%d\n",n,deriv_num,DDS_tol,num_points);
    fprintf(fileID,"error,FW,FWstd,FWsuccess,FWEl,FWElstd,FWElsuccess,AGM,AGMstd,AGMsuccess,AGMEl,AGMElstd,AGMElsuccess\n");
    for i = 1:length(ratio_list)
        fprintf(fileID,"%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",ratio_list(i)*100,...
            mean(t_FW(i,:)'./runtime_DDS*100,"omitnan"),std(t_FW(i,:)'./runtime_DDS*100,"omitnan"),success_FW(i)/num_points*100,...
            mean(t_FW_ele(i,:)'./runtime_DDS*100,"omitnan"),std(t_FW_ele(i,:)'./runtime_DDS*100,"omitnan"),success_FW_ele(i)/num_points*100,...
            mean(t_AGM(i,:)'./runtime_DDS*100,"omitnan"),std(t_AGM(i,:)'./runtime_DDS*100,"omitnan"),success_AGM(i)/num_points*100,...
            mean(t_AGM_ele(i,:)'./runtime_DDS*100,"omitnan"),std(t_AGM_ele(i,:)'./runtime_DDS*100,"omitnan"),success_AGM_ele(i)/num_points*100);
    end
    fprintf(fileID,"\n");


    fclose(fileID);
    
end


end