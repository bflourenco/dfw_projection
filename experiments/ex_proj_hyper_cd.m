%
% Code for the experiments regarding the c_D constant.
%
% Run with n = 10, deriv_num = 1 and with n = 20, deriv_num = 2 in order to 
% generate the files for Table 8 in the paper.
%
clear;
n = 10;
deriv_num = 1;
%Uncomment the next two lines for the experiments in Table 8 - c,d.
%n = 20; 
%deriv_num = 2;

num_points = 30;
DDS_tol = 1e-8;
ratio_list = [1e-1,5e-2,1e-2,5e-3,1e-3,5e-4,1e-4,5e-5,1e-5];
cd_mult = [1,2,4,8,16,100];
cd_mult_size = size(cd_mult,2);
save_to_file = true;
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

% %% measure running time of dual-Frank-Wolfe method
% fprintf('\n start dual-Frank-Wolfe method\n\n')
% % set up
% clearvars poly 
% e = ones(n,1);
% 
% % make poly for hyperbolicity cone constraint
% combi = nchoosek(1:n,n-deriv_num);
% row = size(combi,1);
% poly{1}.type = "matrix";
% poly{1}.mat = zeros(row,n+1);
% for i = 1:row
%     for j = combi(i,:)
%         poly{1}.mat(i,j) = 1;
%     end
% end
% poly{1}.mat(:,n+1) = ones(row,1);
% 
% %%measure runtime
% %prepare data strorages for all datasets d
% x_FW = cell(num_points,cd_mult_size);
% y_FW = cell(num_points,cd_mult_size);
% itr_FW = zeros(num_points,cd_mult_size);
% gap_FW = cell(1,num_points,cd_mult_size);
% feas_FW = cell(1,num_points,cd_mult_size);
% runtime_FW = cell(1,num_points,cd_mult_size);
% obj_vals_FW = cell(1,num_points,cd_mult_size);
% 
% %measure runtime for each d
% for k = 1:num_points
%     fprintf('\n %d-th point\n\n',k)
%     for l = 1:cd_mult_size
%         fprintf('cd multiplier = %g\n',cd_mult(l));
%         clear("opts");
%         opts.y0 = zeros(n,1);
%         %make problem
%         b = -d(:,k);
%         step.rule = "Lipschitz";
%         step.L = 1;
%         opts.step = step;
%         grad = @(x) (x-b);
%         %c = norm(e)*sqrt(norm(e-d(:,k)).^2+norm(e).^2);
%         c = norm(e)*norm(e-d(:,k));
%         if dot(e,opts.y0) > c, c = dot(e,opts.y0); end
%         c = c*cd_mult(l);
%         
%         %opts.zero_eps = 1e-1;
%         %opts.verbose_freq = 10;
%     
%         %stopping criteria
%         opts.stop.maxitr = Inf;
%         opts.stop.feas_check = true;
%     
%         %It is necessary to stop the algorithm when the FW gap is too small or is
%         %zero, since this may lead to numerical problems.
%         %Still, we want to algorithm to keep going until we reach the desired 
%         %precision.
%         %Setting it to 1e-12 ensures that it will rarely (or never) stop 
%         %because of the FW gap stopping criterion.
%         opts.stop.gap_tol = 1e-12;
%     
%         %in order to save time, we stop the experiment if an iterate is 
%         %already within min(ratio_list) of the objective value obtained by DDS.
%         opts.stop.f = @(x) (0.5*norm(x-d(:,k)).^2);
%         opts.stop.obj_thresh = min(obj_vals_DDS{k})*(1+min(ratio_list));
%         opts.stop.max_time = runtime_DDS(k); 
%       
%         [x_FW{k,l},y_FW{k,l},itr_FW(k,l),gap_FW{k,l},feas_FW{k,l},runtime_FW{k,l}] = FW_HP_exp(grad,eye(n),zeros(n,1),poly,e,c,opts);
%         
%         %calculate objective values
%         obj_vals_FW{k,l} = 0.5*vecnorm(x_FW{k,l}-d(:,k)).^2;
%     end
% end

%% measure running time of dual-Frank-Wolfe method using eleSym
fprintf('\n start dual-Frank-Wolfe method using eleSym\n\n')
% remake poly using "symmetric" option
clearvars poly
poly{1}.type = "symmetric";
poly{1}.deg = n-deriv_num;
poly{1}.index = 1:n;

%%measure runtime
%prepare data storages for all datasets d
x_FW_ele = cell(num_points,cd_mult_size);
y_FW_ele = cell(num_points,cd_mult_size);
itr_FW_ele = zeros(num_points,cd_mult_size);
gap_FW_ele = cell(num_points,cd_mult_size);
feas_FW_ele = cell(num_points,cd_mult_size);
runtime_FW_ele = cell(num_points,cd_mult_size);
obj_vals_FW_ele = cell(num_points,cd_mult_size);

%measure runtime for each d
for k = 1:num_points
    fprintf('\n %d-th point\n\n',k)
    for l = 1:cd_mult_size
        fprintf('cd multiplier = %g\n',cd_mult(l));
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
        c = c*cd_mult(l);

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
      
        [x_FW_ele{k,l},y_FW_ele{k,l},itr_FW_ele(k,l),gap_FW_ele{k,l},feas_FW_ele{k,l},runtime_FW_ele{k,l}] = FW_HP_exp(grad,eye(n),zeros(n,1),poly,e,c,opts);
        
        %calculate objective values
        obj_vals_FW_ele{k,l} = 0.5*vecnorm(x_FW_ele{k,l}-d(:,k)).^2;
    end
end



%% calculate mean time to achieve some relative errors
t_FW = NaN(cd_mult_size,length(ratio_list),num_points);
t_FW_ele = NaN(cd_mult_size,length(ratio_list),num_points);
t_FW_it = NaN(cd_mult_size,length(ratio_list),num_points);
t_FW_ele_it = NaN(cd_mult_size,length(ratio_list),num_points);
t_FW_xabs = NaN(cd_mult_size,length(ratio_list),num_points);
t_FW_ele_xabs = NaN(cd_mult_size,length(ratio_list),num_points);
success_FW = zeros(cd_mult_size,length(ratio_list));
success_FW_ele = zeros(cd_mult_size,length(ratio_list));
for l=1:cd_mult_size
    for i = 1:length(ratio_list)
        for k = 1:num_points
%             if ~isempty(runtime_FW{k,l}( (obj_vals_FW{k,l}-obj_vals_DDS{k})<=obj_vals_DDS{k}*ratio_list(i) & feas_FW{k,l} ))
%                 %find the first iteration for which the iterate satisfies
%                 %the given conditions
%                 min_it = find(((obj_vals_FW{k,l}-obj_vals_DDS{k})<=obj_vals_DDS{k}*ratio_list(i)) & feas_FW{k,l},1);
%                 t_FW(l,i,k) = runtime_FW{k,l}(min_it);
%                 t_FW_it(l,i,k) = min_it;
%                 t_FW_xabs(l,i,k) = norm(x_DDS{k}-x_FW{k}(:,min_it),"inf");
%                 %t_FW(i,k) = min(runtime_FW{k}((obj_vals_FW{k}-obj_vals_DDS{k})<=obj_vals_DDS{k}*ratio_list(i) & feas_FW{k}) );
%                 success_FW(l,i) = success_FW(l,i)+1;
%             end
            if ~isempty(runtime_FW_ele{k,l}((obj_vals_FW_ele{k,l}-obj_vals_DDS{k})<=obj_vals_DDS{k}*ratio_list(i) & feas_FW_ele{k,l}))
                min_it = find(((obj_vals_FW_ele{k,l}-obj_vals_DDS{k})<=obj_vals_DDS{k}*ratio_list(i)) & feas_FW_ele{k,l},1);
                t_FW_ele(l,i,k) = runtime_FW_ele{k,l}(min_it);
                t_FW_ele_it(l,i,k) = min_it;
                t_FW_ele_xabs(l,i,k) = norm(x_DDS{k}-x_FW_ele{k,l}(:,min_it),"inf");
                success_FW_ele(l,i) = success_FW_ele(l,i)+1;
            end
        end
    end
end
% display results
fprintf("\n Average time to achieve relative errors with respect to DDS\n\n")
for l=1:cd_mult_size
    fprintf('cd multiplier = %g\n',cd_mult(l));
    fprintf("rel error;\tFW_ele;it;std;rel_t;std;success;\n")    
    for i = 1:length(ratio_list)
%         fprintf("%.2e;\t%.2e;%.2e;%4.2f;%4.2f;%.2e;%.2e;%4.2f;%4.2f;%4.1f;\t%.2e;%.2e;%4.2f;%4.2f;%.2e;%.2e;%4.2f;%4.2f;%4.1f\n",ratio_list(i),...
%             mean(t_FW_xabs(l,i,:),"omitnan"),std(t_FW_xabs(l,i,:),"omitnan"),mean(t_FW_it(l,i,:),"omitnan"),std(t_FW_it(l,i,:),"omitnan"),mean(t_FW(l,i,:),"omitnan"),std(t_FW(l,i,:),"omitnan"),mean(t_FW(l,i,:)'./runtime_DDS*100,"omitnan"),std(t_FW(l,i,:)'./runtime_DDS*100,"omitnan"),success_FW(l,i)/num_points*100,...
%             mean(t_FW_ele_xabs(l,i,:),"omitnan"),std(t_FW_ele_xabs(l,i,:),"omitnan"),mean(t_FW_ele_it(l,i,:),"omitnan"),std(t_FW_ele_it(l,i,:),"omitnan"),mean(t_FW_ele(l,i,:),"omitnan"),std(t_FW_ele(l,i,:),"omitnan"),mean(t_FW_ele(l,i,:)'./runtime_DDS*100,"omitnan"),std(t_FW_ele(l,i,:)'./runtime_DDS*100,"omitnan"),success_FW_ele(l,i)/num_points*100);          
         fprintf("%.2e;\t%4.2f;%4.2f;%4.2f;%4.2f;%4.1f\n",ratio_list(i),...
            mean(t_FW_ele_it(l,i,:),"omitnan"),std(t_FW_ele_it(l,i,:),"omitnan"),mean(reshape(t_FW_ele(l,i,:),[1,num_points])'./runtime_DDS*100,"omitnan"),std(reshape(t_FW_ele(l,i,:),[1,num_points])'./runtime_DDS*100,"omitnan"),success_FW_ele(l,i)/num_points*100);          
 
    
    end
    fprintf('\n');
end

if save_to_file
    filename = "proj_hyper_cd_"+sprintf("n%d_",n)+sprintf("d%d",deriv_num)+sprintf("_%d",num_points);
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
    fprintf(fileID,"error,");
    for l=1:(cd_mult_size-1)
        prefix = "cd"+string(l);
        fprintf(fileID,prefix+"FWElit,"+prefix+"FWElitstd,"+prefix+"FWElabs,"+prefix+"FWElabsstd,"+prefix+"FWEl,"+prefix+"FWElstd,"+prefix+"FWElsuccess,");
    end
    prefix = "cd"+string(cd_mult_size);
    fprintf(fileID,prefix+"FWElit,"+prefix+"FWElitstd,"+prefix+"FWElabs,"+prefix+"FWElabsstd,"+prefix+"FWEl,"+prefix+"FWElstd,"+prefix+"FWElsuccess\n");

    for i = 1:length(ratio_list)
        fprintf(fileID,"%g,",ratio_list(i)*100);
        for l=1:(cd_mult_size-1)
            fprintf(fileID,"%g,%g,%.2e,%.2e,%g,%g,%g,",...
              mean(t_FW_ele_it(l,i,:),"omitnan"),std(t_FW_ele_it(l,i,:),"omitnan"),mean(t_FW_ele(l,i,:),"omitnan"),std(t_FW_ele(l,i,:),"omitnan"),mean(reshape(t_FW_ele(l,i,:),[1,num_points])'./runtime_DDS*100,"omitnan"),std(reshape(t_FW_ele(l,i,:),[1,num_points])'./runtime_DDS*100,"omitnan"),success_FW_ele(l,i)/num_points*100);          
        end
        l = cd_mult_size;
        fprintf(fileID,"%g,%g,%.2e,%.2e,%g,%g,%g\n",...
              mean(t_FW_ele_it(l,i,:),"omitnan"),std(t_FW_ele_it(l,i,:),"omitnan"),mean(t_FW_ele(l,i,:),"omitnan"),std(t_FW_ele(l,i,:),"omitnan"),mean(reshape(t_FW_ele(l,i,:),[1,num_points])'./runtime_DDS*100,"omitnan"),std(reshape(t_FW_ele(l,i,:),[1,num_points])'./runtime_DDS*100,"omitnan"),success_FW_ele(l,i)/num_points*100);          
   
    end
    fprintf(fileID,"\n");


    fclose(fileID);
    
end


