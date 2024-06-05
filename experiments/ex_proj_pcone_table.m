function[] = ex_proj_pcone_table(n,p,num_points,IPM_tol,ratio_list,skip_DDS,save_to_file)


% make the minimum eigenvalue and normal vector oracle
mineig = @(x,Te) (x(1)-norm(x(2:end),p));
normal = @(x) ([1;-sign(x(2:end)).*(abs(x(2:end)/x(1)).^(p-1))]);
e = [1;zeros(n,1)];

% Generate points to project
rng(1); %fix seed
d = zeros(n+1,num_points);
for k = 1:num_points
    dk = randn(n+1,1);
    min_eig_dk = mineig(dk,e);
    %Discard test points that are too close to the cone.
    while (min_eig_dk > -1e-4 )
        dk = randn(n,1);
        min_eig_dk = mineig(dk,e);
        fprintf("Discarded a test point because it is too close to being feasible\n")
    end
    fprintf("Min eigenvalue of the %dth point: %g \n",k,min_eig_dk)
    d(:,k) = dk;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% measure running time of DDS
fprintf('\n start DDS\n\n')

% set up
runtime_DDS = zeros(num_points,1);%cell(1,num_points);
x_DDS = cell(1,num_points);
obj_vals_DDS = cell(1,num_points);
itr_DDS = cell(1,num_points);


% measure running time
for k = 1:num_points
    if skip_DDS == false
            fprintf('\n %d-th point\n\n',k)
            %%%%% start making problem %%%%%
            % adding SOCP constraint
            cons{1,1} = 'SOCP';
            cons{1,2} = n+1;
            A{1,1} = [1,zeros(1,2*n+1);
                      zeros(n+1,1),eye(n+1),zeros(n+1,n)];
            b{1,1} = [0;-d(:,k)];
         
            % adding LP constraint
            cons{2,1}= 'LP';
            cons{2,2}= 2;
            A{2,1} = [0,1,zeros(1,2*n);
                      0,1,zeros(1,n),-ones(1,n)];
            b{2,1} = [0;0];
         
            % adding p-norm constraint
            cons{3,1}='GPC';
            A{3,1}=[];
            for i = 1:n
                cons{3,2}{i,1} = [2,1];
                cons{3,2}{i,2} = [1-(1/p);1/p];
                A{3,1} = [A{3,1};
                          0,1,zeros(1,2*n);
                          zeros(1,n+i+1),1,zeros(1,n-i);
                          zeros(1,1+i),1,zeros(1,2*n-i)];
            end
            b{3,1} = zeros(3*n,1);
         
            % set cost vector c
            c=[1;zeros(2*n+1,1)];
            %%%%% finish making problem %%%%%
        
            OPTIONS.tol = IPM_tol;
            info.status = 1;
            [x,~,info]=DDS(c,A,b,cons,OPTIONS);
            x_DDS{k} = x(2:n+2);
            runtime_DDS(k) = info.time;
            itr_DDS{k} = info.iter;
            obj_vals_DDS{k} = 0.5*vecnorm(x_DDS{k}-d(:,k)).^2;
    else
        obj_vals_DDS{k} = nan;
        runtime_DDS(k) = nan;
    end
    

end

%% measure running time of Mosek
fprintf("\n start Mosek\n\n")
%set up
clearvars prob info param

% Prepare storages for results
runtime_Mosek = zeros(num_points,1); %cell(1,num_points);
x_Mosek = cell(1,num_points);
itr_Mosek = cell(1,num_points);
obj_vals_Mosek = cell(1,num_points);

for k = 1:num_points
    %%%%% start making problem %%%%%
    [~, res] = mosekopt('symbcon');

    % Specify the non-conic part of the problem.
    % Variables number 1,2~n+2,n+3~2*n+2 correspond to t,x,y
    prob.c = [1,zeros(1,2*n+1)];
    prob.a = [0,1,zeros(1,n),-ones(1,n)];
    prob.blc = [0];
    prob.buc = [Inf];
    prob.blx = [-Inf, 0,-Inf(1,2*n)];
    prob.bux = [Inf(1,2*n+2)];

    % Specify the cones as affine conic constraints.
    % adding the Second-order cone constraint, which is dependent on d(:,k)
    prob.f = sparse([eye(n+2),zeros(n+2,n)]);
    prob.g = [0;-d(:,k)];
    prob.accs = [res.symbcon.MSK_DOMAIN_QUADRATIC_CONE n+2];

    % adding the power cone constraints
    for i = 1:n
        prob.f = [prob.f;
                  sparse([0,1,zeros(1,2*n);
                  zeros(1,1+n+i),1,zeros(1,n-i);
                  zeros(1,1+i),1,zeros(1,2*n-i)])];
        prob.g = [prob.g;zeros(3,1)];
        prob.accs = [prob.accs, res.symbcon.MSK_DOMAIN_PRIMAL_POWER_CONE 3 2 1-(1/p) (1/p)];
    end
    %%%%% finish making problem %%%%%


    param.MSK_DPAR_INTPNT_CO_TOL_MU_RED = IPM_tol;
    [~,res]=mosekopt('minimize info',prob, param);
    x_Mosek{k}(:,end+1) = res.sol.itr.xx(2:2+n);
    runtime_Mosek(k) = res.info.MSK_DINF_OPTIMIZER_TIME;
    itr_Mosek{k}(end+1) = res.info.MSK_IINF_INTPNT_ITER;
 
    obj_vals_Mosek{k} = 0.5*vecnorm(x_Mosek{k}-d(:,k)).^2;
end

%% measure running time of dual-Frank-Wolfe method
fprintf('\n start dual-Frank-Wolfe method\n\n')

% set up
clearvars opts


% prepare data strorage for results
itr_FW = zeros(1,num_points);
gap_FW = cell(1,num_points);
feas_FW = cell(1,num_points);
runtime_FW = cell(1,num_points);
obj_vals_FW = cell(1,num_points);
x_FW = cell(1,num_points);

%measure runtime for each d
for k = 1:num_points
    fprintf('\n %d-th point\n\n',k)
    clear("opts");
    %make problem
    step.L = 1;
    step.rule = "Lipschitz";
    opts.y0 = zeros(n+1,1);
    opts.step = step;
    obj_func = @(x) 0.5*(norm(x-d(:,k))).^2;
    grad = @(x) (x+d(:,k));
    c = norm(e)*norm(e-d(:,k));
    if dot(e,opts.y0) > c, c = dot(e,opts.y0); end
    
    %stopping criteria
    opts.stop.maxitr = Inf;
    opts.stop.gap_tol = 1e-16;
    opts.stop.feas_check = true;
    opts.stop.f = obj_func;
    opts.stop.max_time = runtime_Mosek(k); %runtime_DDS(k); 
    
    [x_FW{k},~,itr_FW(k),gap_FW{k},feas_FW{k},runtime_FW{k}] = FW_GCP_exp(grad,eye(n+1),zeros(n+1,1),e,c,normal,mineig,opts);
    obj_vals_FW{k} = 0.5*vecnorm(x_FW{k}-d(:,k)).^2;
    %[obj_vals_FW{k},itr_FW(k),gap_FW{k},feas_FW{k},runtime_FW{k}] = FW_GCP_LessMemory(obj_func,grad,eye(n+1),zeros(n+1,1),e,c,normal,mineig,opts);
end

%% calculate mean time to achieve some relative errors

t_FW_Mosek = NaN(length(ratio_list),num_points);
t_FW_DDS = NaN(length(ratio_list),num_points);
success_FW_Mosek = zeros(length(ratio_list),1);
success_FW_DDS = zeros(length(ratio_list),1);
for k = 1:num_points
    for i = 1:length(ratio_list)
        if ~isempty(runtime_FW{k}(feas_FW{k}&((obj_vals_FW{k}-obj_vals_Mosek{k})<=obj_vals_Mosek{k}*ratio_list(i))))
            success_FW_Mosek(i) = success_FW_Mosek(i)+1;
            t_FW_Mosek(i,k) = min(runtime_FW{k}((obj_vals_FW{k}-obj_vals_Mosek{k})<=obj_vals_Mosek{k}*ratio_list(i)&feas_FW{k}));
        end
        if ~isempty(runtime_FW{k}(feas_FW{k}&((obj_vals_FW{k}-obj_vals_DDS{k})<=obj_vals_DDS{k}*ratio_list(i))))
            success_FW_DDS(i) = success_FW_DDS(i)+1;
            t_FW_DDS(i,k) = min(runtime_FW{k}((obj_vals_FW{k}-obj_vals_DDS{k})<=obj_vals_DDS{k}*ratio_list(i)&feas_FW{k}));
        end
    end
end

fprintf("\n Average time to achieve relative errors with respect to DDS\n\n")
fprintf("relative error;\tFW (Mosek);std;success;\t\tFW (DDS);std;success;\n")
for i = 1:length(ratio_list)
    fprintf("%.2e;\t%4.2f;%4.2f;%4.1f;\t%4.2f;%4.2f;%4.1f\n",ratio_list(i),...
        mean(t_FW_Mosek(i,:)'./runtime_Mosek*100,"omitnan"),std(t_FW_Mosek(i,:)'./runtime_Mosek*100,"omitnan"),success_FW_Mosek(i)/num_points*100,...
        mean(t_FW_DDS(i,:)'./runtime_DDS*100,"omitnan"),std(t_FW_DDS(i,:)'./runtime_DDS*100,"omitnan"),success_FW_DDS(i)/num_points*100)
end

if save_to_file
    filename = "proj_pcone_"+sprintf("n%d_",n)+replace(sprintf("p%.1f",p),".","_")+sprintf("_%d",num_points);
    if (IPM_tol < 1e-3)
        filename = filename + sprintf("_tol_high");
    else
        filename = filename + sprintf("_tol_low");
    end
    datetime.setDefaultFormats('default','yyyy_MM_dd_hh_mm_ss');
    filename = filename + "_" + string(datetime("now"));
    save(filename + ".mat");
    datetime.setDefaultFormats('reset');
    
    fileID = fopen(filename +"stats" +".csv",'w');
    %fprintf(fileID,"n;p;tol;points\n");
    %fprintf(fileID,"%d;%g;%.1e;%d\n",n,p,IPM_tol,num_points);
    if (skip_DDS==true)
            fprintf(fileID,"error,FWM,FWMstd,FWMSuccess\n");
            for i = 1:length(ratio_list)
                fprintf(fileID,"%g,%g,%g,%g\n",ratio_list(i)*100,...
                    mean(t_FW_Mosek(i,:)'./runtime_Mosek*100,"omitnan"),std(t_FW_Mosek(i,:)'./runtime_Mosek*100,"omitnan"),success_FW_Mosek(i)/num_points*100);
            end
    else
            fprintf(fileID,"error,FWM,FWMstd,FWMSuccess,FWD,FWDstd,FWDsuccess\n");
            for i = 1:length(ratio_list)
                fprintf(fileID,"%g,%g,%g,%g,%g,%g,%g\n",ratio_list(i)*100,...
                    mean(t_FW_Mosek(i,:)'./runtime_Mosek*100,"omitnan"),std(t_FW_Mosek(i,:)'./runtime_Mosek*100,"omitnan"),success_FW_Mosek(i)/num_points*100,...
                     mean(t_FW_DDS(i,:)'./runtime_DDS*100,"omitnan"),std(t_FW_DDS(i,:)'./runtime_DDS*100,"omitnan"),success_FW_DDS(i)/num_points*100);
            end
    end

    %fprintf(fileID,"\n");


    fclose(fileID);    
end
end
