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
function[x,y,itr,gap,feas,time] = FW_GCP_exp(grad,T,b,e,c_D,normal_oracle,mineig_oracle,opts)
%This function solves the following problem by Frank-Wolfe method.
%min  f(x)
%s.t. Tx+u in Lambda(p,e)
%where  f : (R^n) to (R or {infty}) is a mu-strongly and closed proper convex function
%       T : (R^n) to  (R^m) is a linear map
%       b is included in (R^m)
%       Lambda(p,e) (included in R^m) is a pointed hyperbolicity cone associated with hyperbolic polynomial poly
%       e (in R^n) is a vector satisfying e in ri(Lambda(p,e))

%Input[datatype]
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
%              gap_tol - tolerance of FW gap [real number] (1e-2)
%              maxitr - max number of iteration [positive integer] (5000)
%              max_time - max time / unit is second [real number] (Inf)
%              obj_thresh - tolerance of objective value [real number] (-Inf)
%                        opts.stop.f must be given to specify this variable.
%              f - objective function [function handle]

%Output [datatype]
%x - points in primal space   x(:,i) corresponds to x_i in algorithm [m by iter matrix]
%s - points in dual space  s(:,i) corresponds to s_i in algorithm [m by iter matrix]
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
if isfield(opts,'step'), step = opts.step; else step.rule = 'deminishing'; end
if isfield(opts,'zero_eps'), zero_eps = opts.zero_eps; else zero_eps = 1e-8; end
if isfield(opts,'verbose_freq')
    verbose_freq = opts.verbose_freq;
    fprintf(" Iteration          Time              FW_gap     \n")
    fprintf("-------------------------------------------------\n")
else verbose_freq= Inf; end
% stopping criteria
if isfield(opts,'stop') == false, opts.stop = {}; end
if isfield(opts.stop,'gap_tol'), gap_tol = opts.stop.gap_tol; else gap_tol = 1e-2; end
if isfield(opts.stop,'maxitr'), maxitr = opts.stop.maxitr; else maxitr = 5000; end
if isfield(opts.stop,'max_time'), max_time = opts.stop.max_time; else max_time = Inf; end
if isfield(opts.stop,'f'), f = opts.stop.f; else f = @(x) (Inf); end
if isfield(opts.stop,'obj_thresh')
    if ~isfield(opts.stop,'f')
        error('To specify opts.stop.obj_val, set opts.stop.f.');
    end
    obj_thresh = opts.stop.thresh;
else
    obj_thresh = -Inf;
end


%%set up
TT = T.';
gap = zeros(1,1); %gap is the list of Frank-Wolfe gap
y(:,1) = y0; %s is the list of the dual variables
x = zeros(n,1); %s is the list of the primal variables
time = zeros(1,1); %time is the list of the time
feas = true(1,1); %feas is the list of the feasibility of x

%%main algorithm
itr = 1;
while(1)
    x(:,itr) = grad(TT*y(:,itr));
    g = T * x(:,itr) + b;
    
    %calculate the minimum eigenvalue of g
    mineig_g = mineig_oracle(g,e);
    
    %solve subproblem and calculate descent direction
    if mineig_g > -zero_eps % If g is included in Lambda(p,e)
        feas(itr) = true;
        descent = -y(:,itr);
    else
        feas(itr) = false;
        z = g - mineig_g*e;
        normal_vec = normal_oracle(z);
        descent = c_D*normal_vec/dot(e,normal_vec) - y(:,itr);
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
        f_k = f(x(:,itr)); 
        if ( gap(itr) < gap_tol || f_k < obj_thresh), break, end
    end
    
    %update
    alpha = stepsize(step,y(:,itr),g,descent,itr); %alpha is the stepsize determined by stepsize_rule
    if alpha > 1, alpha = 1; end
    y(:,itr+1) = y(:,itr) + alpha*descent;
    itr = itr + 1;
end
x = x(:,1:itr);
y = y(:,1:itr);
gap = gap(1:itr);
feas = feas(1:itr);
end