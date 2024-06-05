% Copyright (C) 2024 Takayuki Nagano, Bruno F. Lourenço, Akiko Takeda
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
function[obj_vals,itr,gap,feas,time] = FW_GCP_LessMemory(f,grad,T,b,e,c_D,normal_oracle,mineig_oracle,opts)
%This function solves the following problem by Frank-Wolfe method.
%min  f(x)
%s.t. Tx+b in Lambda(p,e)
%where  f : (R^n) to (R or {infty}) is a mu-strongly and closed proper convex function
%       T : (R^n) to  (R^m) is a linear map
%       b is included in (R^m)
%       Lambda(p,e) (included in R^m) is a pointed hyperbolicity cone associated with hyperbolic polynomial poly
%       e (in R^n) is a vector satisfying T*e in ri(Lambda(p,e))
%This function returns only {f(x_k)} and G(y_k). It does not return {x_k}
%and {s_k}. If the dimention is large, {x_k} and {s_k} requires huge memory
%and FW_GCP.m would not run. In such cases, this function is recommended.

%Input[datatype]
%f - the objective function [function handle]
%grad - the gradient of f^* [function handle from (n by 1 vector) to (n by 1 vector)]
%T - [m by n matrix]
%b - [m by 1 vector]
%e - [n by 1 vector]
%c_D - real number satisfying Assumption 1 [real number]
%normal - function handle which return a normal vector of \Lambda at x [function handle from (n by 1 vector) to (n by 1 vector)]
%mineig_oracle - mineig_oracle is the function handle which calculates the minimum eigenvalues of x.
%                [function handle from (n by 1 vector) to (n by 1 vector)]
%                The arguments of mineig_oracle is (x,e).
%opts - contains the following structure variables[datatype], default values are in ():
%       y0 - starting point [m by 1 vector] (0)
%       feas_check - [logical] (true) If feas_check is true, the stopping
%                    criteria are judged only when x_i is included in feasible region.
%       step - specifying step size rule[structure variables](deminishing step size rule).
%              Detailed information is in stepsize.m
%       zero_eps - margin used for judging equalities and inequalities. [real number] (1e-8)
%       verbose_freq - [positive integer] (Inf)
%                      If verbose_freq is specified, some informations are
%                      displayed every opts.verbose_freq iterations.
%       eig_oracle{1} - eig_oracle{1}(x,e) is the function which calculates the eigenvalues of x w.r.t. poly{1} and e.
%                       If you give this oracle, this program calculates eigenvalues using analytic form. 
%                       Otherwise, this program calculates them by finding the roots of characteristic polynomial.
%                       You can also input eig_oracle_2,3,... for poly{2},poly{3},...
%       stop - parameters for setting stopping criteria [structure variables]
%              If one of the criteria described below was achieved, this
%              program stops.
%              maxitr - max number of iteration [positive integer] (5000)
%              max_time - max time / unit is second [real number] (Inf)
%              gap_tol - tolerance of FW gap [real number] (1e-2)
%              obj_thresh - tolerance of objective value [real number] (-Inf)
%              dual_gap - tolerance of duality gap [real number] (Inf)
%                         opts.stop.f_conj must be given to specify this variable.
%              KKT_tol - tolerance of the sum of duality gap and min{0,lambda_min(Tx_k+b)} [real number] (Inf)
%                        opts.stop.f_conj must be given to specify this variable.
%              f_conj - conjugate function of objective function [function handle]

%Output[datatype]
%obj_vals - the objective values f(x_k) [1by iter matrix]
%itr - the number of iteration actually done in the program [positive integer]
%gap - Frank-Wolfe gap gap(i) corresponds to Frank-Wolfe gap at s_i [m by iter matrix]
%feas - feas(i) is 1 if x_i is feasible and 0 otherwise [1 by iter matrix]
%time - time ellapsed after the i-th iteration [1 by iter matrix]

%%starting point of time measuring
tic;

%%preprocessing options
[m,n] = size(T);
if isfield(opts,'y0'), y0 = opts.y0; else y0 = zeros(m,1); end 
if isfield(opts,'feas_check'), feas_check = opts.feas_check; else feas_check = true; end
if isfield(opts,'step'), step = opts.step; else step.rule = 'diminishing'; end
if isfield(opts,'zero_eps'), zero_eps = opts.zero_eps; else zero_eps = 1e-8; end
if isfield(opts,'verbose_freq')
    verbose_freq = opts.verbose_freq;
    fprintf(" Iteration          Time              FW_gap     \n")
    fprintf("-------------------------------------------------\n")
else verbose_freq= Inf; 
end
% stopping criteria
if isfield(opts,'stop') == false, opts.stop = {}; end
if isfield(opts.stop,'gap_tol'), gap_tol = opts.stop.gap_tol; else gap_tol = 1e-2; end
if isfield(opts.stop,'maxitr'), maxitr = opts.stop.maxitr; else maxitr = 5000; end
if isfield(opts.stop,'max_time'), max_time = opts.stop.max_time; else max_time = Inf; end
if isfield(opts.stop,'f_conj'), f_conj = opts.stop.f_conj; else f_conj = @(x) (Inf); end
if isfield(opts.stop,'obj_thresh')
    if ~isfield(opts.stop,'f')
        error('To specify opts.stop.obj_val, set opts.stop.f.');
    end
    obj_thresh = opts.stop.obj_thresh;
else
    obj_thresh = -Inf;
end
if isfield(opts.stop,'dual_gap')
    if ~(isfield(opts.stop,'f')&&isfield(opts.stop,'f_conj'))
        error('To specify opts.stop.obj_val, set opts.stop.f and opts.stop.f_conj.');
    end
    dual_gap = opts.stop.dual_gap;
else
    dual_gap = -Inf;
end
if isfield(opts.stop,'KKT')
    if ~(isfield(opts.stop,'f')&&isfield(opts.stop,'f_conj'))
        error('To specify opts.stop.obj_val, set opts.stop.f and opts.stop.f_conj');
    end
    KKT_tol = opts.stop.KKT_tol;
else
    KKT_tol = -Inf;
end

%%set up
TT = T.';
gap = zeros(1,1); %gap is the list of Frank-Wolfe gap
y = y0; %s is the list of the dual variables
time = zeros(1,1); %time is the list of the time
obj_vals = zeros(1,1); %obj_vals is the list of the objective value
feas = true(1,1); %feas is the list of the feasibility of x

%%main algorithm
itr = 1;
while(1)
    x = grad(TT*y);
    obj_vals(itr) = f(x);
    g = T * x + b;
    
    %calculate the minimum eigenvalue of g
    mineig_g = mineig_oracle(g,e);
    
    %solve subproblem and calculate descent direction
    if mineig_g > -zero_eps % If g is included in Lambda(p,e)
        feas(itr) = true;
        descent = -y;
    else
        feas(itr) = false;
        z = g - mineig_g*e;
        normal_vec = normal_oracle(z);
        descent = c_D*normal_vec/dot(e,normal_vec) - y;
    end
    
    %calculate the Frank-Wolfe gap and measure time
    gap(itr) = -dot(g,descent);
    time(itr) = toc;

    %display intermediary results
    if rem(itr,verbose_freq)==0
        str_time = sprintf("%.5f",time(itr));
        str_gap = sprintf("%.5e",gap(itr));
        str_itr = sprintf("%d",itr);
        fprintf(blanks(10-strlength(str_itr))+str_itr+'   |  ' ...
        +blanks(11-strlength(str_time))+str_time+'   |   '+blanks(13-strlength(str_gap))+str_gap+'\n')
    end    
    
    %If one of the stopping criteria is satisfied, stop
    if itr>=maxitr, fprintf("Reached max iteration.\n"), break, end
    if time(end) > max_time, fprintf("Reached max time.\n"), break, end
    if feas(itr) || ~feas_check
        f_k = f(x); f_conj_k = f_conj(y);
        if ( gap(itr) < gap_tol || f_k < obj_thresh || (f_k+f_conj_k) < dual_gap || (f_k+f_conj_k + mineig_g) < KKT_tol ), break, end
    end
    
    %update
    alpha = stepsize(step,y,g,descent,itr); %alpha is the stepsize determined by stepsize_rule
    if alpha > 1, alpha = 1; end
    y = y + alpha*descent;
    itr = itr + 1;
end
gap = gap(1:itr);
feas = feas(1:itr);
end