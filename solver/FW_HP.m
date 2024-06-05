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
function[x,y,status] = FW_HP(grad_f_conj,T,b,poly,e,c_D,opts)
%This function solves the following problem by the dual-Frank-Wolfe method.
%min  f(x)
%s.t. Tx+b in Lambda(p,e)
%where  f : (R^n) to (R or {infty}) is a mu-strongly and closed proper convex function
%       T : (R^n) to  (R^m) is a linear map 
%       b belongs to (R^m)
%       Lambda(p,e) (included in R^m) is a pointed hyperbolicity cone associated with hyperbolic polynomial 'poly'
%       e (in R^m) is a vector in ri(Lambda(p,e))
%

%Input[datatype]
%grad_f_conj - the gradient of f^* [function handle from (n by 1 vector) to (n by 1 vector)]
%T - [m by n matrix]
%u - [m by 1 vector]
%poly - the set of polynomials from R^m to R [structure variables]
%       poly{1}*...*poly{end} represents the whole polynomial.
%       poly{i} is represented by the matrix or oracle, and special option is prepared for elementary symmetric polynomial.
%       Each representation requires following items.
%       <matrix representation>
%       To use matrix representation, set poly{i}.type "matrix".
%       poly{i}.mat - [(the number of terms) times n matrix] Each row represents a term,
%                     and elements in the row correspond to the exponents of variables.
%       <oracle representation>
%       To use oracle representation, set poly{i}.type "oracle".
%       poly{i}.p - [function handle] oracle for evaluating poly{i} at x
%       poly{i}.deg - [integer] degree of poly{i}
%       poly{i}.grad -[function handle] oracle for evaluating the gradient of poly{i} at x
%       <elementary symmetric polynomial>
%       If (x_i1,...x_ik) is in the hyperbolicity cone associated with the k-th symmetric polynomial,
%       set poly{i}.type = "symmetric", poly{i}.index = [i1,...,ik], and poly{i}.deg = k. If you use this option,
%       (e_i1,...,e_ik) must be (1,...,1).
%       
%       Optionally you can give the eigenvalues oracle for poly{i} in
%       poly{i}.eig. 
%e - [m by 1 vector]
%c_D - real number satisfying Assumpution 1 [real number]
%opts - contains the following structure variables[datatype], default values are in ():
%       y0 - starting point [m by 1 vector] (0)
%       step - specifying step size rule[structure variables](diminishing step size rule).
%              Detailed information is in stepsize.m
%       zero_eps - margin used for judging equalities and inequalities. [real number] (1e-8)
%       verbose_freq - [positive integer] (Inf)
%                      If verbose_freq is specified, some informations are
%                      displayed every opts.verbose_freq iterations.
%       stop - parameters for setting stopping criteria [structure variables]
%              If one of the criteria described below was achieved, this
%              program stops.
%              gap_tol - tolerance of FW gap [real number] (1e-2)
%              maxitr - max number of iteration [positive integer] (Inf)
%              max_time - max time / unit is second [real number] (10)
%              obj_thresh - tolerance of objective value [real number] (-Inf)
%                        opts.stop.f must be given to specify this variable.
%              f - objective function [function handle]

%Output[datatype]
%x - A primal solution corresponding to the best dual solution found. [real vector]
%y - Best dual solution found: among the iterates that are primal feasible within opt.zero_eps,
% return the dual solution with minimal FW gap. If no such iterate exists, y is
% the dual solution corresponding to a primal solution found with minimal feasibility violation.
%status - contains a structure with the following values:
%
%       status - 1: Obtained a feasible solution with FW gap <= stop.gap_tol
%                2: Obtained a feasible solution but FW gap stop criterion
%                was not satisfied.
%                3: No primal feasible solution found
%       time: total running time.
%       itr: total number of iterations.
%       FW_gap: Best FW_gap found. It is Inf if no primal feasible solution
%       was found
%
%

%% start of running time measurement
tic;

%% preprocessing

% options
[m,n] = size(T);
if isfield(opts,'y0'), y0 = opts.y0; else y0 = zeros(m,1); end 
if isfield(opts,'step'), step = opts.step; else step.rule = 'diminishing'; end
if isfield(opts,'zero_eps'), zero_eps = opts.zero_eps; else zero_eps = 1e-8; end
if isfield(opts,'verbose_freq')
    verbose_freq = opts.verbose_freq;
    fprintf(" Iteration          Time              FW_gap         min_eig    \n")
    fprintf("----------------------------------------------------------------\n")
else verbose_freq= Inf; end
% stopping criteria
if isfield(opts,'stop') == false, opts.stop = {}; end
if isfield(opts.stop,'gap_tol'), gap_tol = opts.stop.gap_tol; else gap_tol = 1e-2; end
if isfield(opts.stop,'maxitr'), maxitr = opts.stop.maxitr; else maxitr = Inf; end
if isfield(opts.stop,'max_time'), max_time = opts.stop.max_time; else max_time = 10; end

% preprocessing the polynomial
whole_deg = 0;
for i = 1:length(poly)
    if poly{i}.type == "symmetric"
        poly{i}.p = @(x) (eleSym(x(poly{i}.index),poly{i}.deg));
        poly{i}.grad = @(x) (grad_eleSym(x,poly{i}.index,poly{i}.deg));

    elseif poly{i}.type == "matrix"
        poly{i}.p = @(x) (sum(prod(cat(1,x.^((poly{i}.mat(:,1:end-1)).'),(poly{i}.mat(:,end)).'))));
        poly{i}.deg = sum(poly{i}.mat(1,1:end-1));
        grad_poly = cell(1,m);
        for term = 1:size(poly{i}.mat,1)
            for j = 1:m
                if poly{i}.mat(term,j) ~= 0
                    tmp = poly{i}.mat(term,:); tmp(m+1) = tmp(j)*tmp(m+1); tmp(j) = tmp(j)-1;
                    grad_poly{j}(end+1,:) = tmp;
                end
            end
        end
        poly{i}.grad = @(x) (poly_vec_val(x,grad_poly));
    end
    whole_deg = whole_deg + poly{i}.deg;  %whole_deg is the degree of polynomial
end

% set up
TT = T.';
gap = Inf; %Frank-Wolfe gap
y = y0;
best_gap = Inf;
best_eig = -Inf;
best_x = zeros(n,1);
best_y = zeros(m,1);
status.status = 3;

%% main algorithm
itr = 1;
while(1)
    x = grad_f_conj(TT*y);
    g = T * x + b;
    
    %calculate eigenvalues of g
    eig_g = eigH(g,e,poly);
    min_eig_g = min(eig_g);
    feas = false;
    
    %solve subproblem and calculate descent direction
    if min_eig_g > -zero_eps % If g is included in Lambda(p,e)
        feas = true;
        descent = -y;
    else
        z = g - min_eig_g*e;
        eig_z = eig_g - min_eig_g; % eig_z is the list of eigenvalues of z
        mult = nnz(abs(eig_z) < zero_eps); %mult is the multiplicity of z
        normal_vec = real(grad_deriv_poly(z,mult-1,poly,whole_deg,e));
        descent = c_D*normal_vec/dot(e,normal_vec) - y;
    end
    
    %calculate the Frank-Wolfe gap 
    gap = -dot(g,descent);
    
    %If one of the stopping criteria is satisfied, stop
    if feas 
        if (gap <= gap_tol) 
            status.status = 1;
            best_y = y;
            best_x = x;
            best_gap = gap;
            break, 
        elseif (gap < best_gap)
            status.status = 2;
            best_y = y;
            best_x = x;
            best_gap = gap;
        end
    elseif status.status == 3 && min_eig_g > best_eig 
            best_eig = min_eig_g;
            best_y = y;
            best_x = x;
    end

    if itr >= maxitr
        fprintf("Reached max iteration.\n");
        break;
    end
    if toc >= max_time
        fprintf("Reached max time.\n");
        break;
    end


    %display intermediary results if verbose_freq is set.
    if rem(itr,verbose_freq)==0
        str_time = sprintf("%.5f",toc);
        str_gap = sprintf("%.5e",gap);
        str_itr = sprintf("%d",itr);
        str_min_eig = sprintf("%.5e",min(eig_g));
        fprintf(blanks(10-strlength(str_itr))+str_itr+'   |  ' ...
        +blanks(11-strlength(str_time))+str_time+'   |   '+blanks(13-strlength(str_gap))+str_gap ...
        +'   |   '+blanks(13-strlength(str_min_eig))+str_min_eig+ '\n')
    end
    
    %update
    alpha = stepsize(step,y,g,descent,itr); % alpha is the stepsize determined by stepsize_rule
    if alpha > 1, alpha = 1; end
    y = y + alpha*descent;
    itr = itr + 1;
end
x = best_x;
y = best_y;
status.FW_gap = best_gap;
status.itr = itr;
status.time = toc;
end