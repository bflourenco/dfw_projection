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
function[points,w1,w2,itr,num_grad_cal,time] = AGM_HP(c,A,v0,poly,e,L,opts)
%This function solves the following problem by Renegar's method.
%min  <c,x>
%s.t. Ax = b
%     x in Lambda(p,e)
%where  c : (in R^n) is a cost vector
%       A : (R^n) to  (R^m) is a linear map 
%       Lambda(p,e) (included in R^n) is a pointed hyperbolicity cone associated with hyperbolic polynomial poly

%Input[datatype]
%c - [n by 1 vector]
%A - [m by n matrix]
%v0 - Initial point which must satisfy A*v0 = b and <c,v0> < <c,e>
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
%e - Direction vector of p which must be included in int(C_p) and satisfy Ae = b [n by 1 vector]
%L - Estimated Lipschitz constant [real value].
%opts - contains the following structure variables[datatype], default values are in ():
%       eps - tolerance [real number] (1e-7)
%       eig_oracle - eig_oracle(x,e) is the function which calculates the eigenvalues of x along direction vector e.
%                    If you give this oracle, this program calculates eigenvalues using analytic form. 
%                    Otherwise, this program calculates them by finding the roots of characteristic polynomial.
%       This program stops when one of the following stopping criteria is satisfied.
%       max_num - maimum number of times of gradient calculation [positive integer] (Inf)
%                 When gradiens are calculated more than this number, this
%                 program stops at the last of the iteration.
%       max_itr - max number of iteration [positive integer] (10)
%       max_time - max time / unit is second (Inf)
%       f_threshold - threshold of objective value [structure variables]
%                     This must contains f [function handle] (f(x) = 0) and
%                     threshold [real number] (-Inf)
%       zero_eps - margin used for judging equalities and inequalities. [real number] (1e-8)


%Output[datatype]
%points - points in the original space [n by iter matrix]
%w1 - points generated by gradient descent method for f_mu where mu = 1/12log(m)
%    w(:,i) corresponds to w_i in algorithm [n by iter matrix]
%w2 - output of gradient descent method for mu = tol/6log(m)]
%itr - the number of iteration actually done in the program [positive integer]
%time - time(i) is the time to finish i-th iteration [1 by iter matrix]


%%starting point of time measuring
tic;

%%options
if isfield(opts,'eps'), eps = opts.eps; else eps = 1e-7; end
if isfield(opts,'max_num'), max_num = opts.max_num; else max_num = Inf; end
if isfield(opts,'max_itr'), max_itr = opts.max_itr; else max_itr = 10; end
if isfield(opts,'max_time'), max_time = opts.max_time; else max_time = Inf; end
if isfield(opts,'f_threshold'), f = opts.f_threshold.f; threshold = opts.f_threshold.threshold;
else f = @(x) (0); threshold = -Inf; end
if isfield(opts,'zero_eps'), zero_eps = opts.zero_eps; else zero_eps = 1e-8; end
n = length(c);

% preprocessing the polynomial
whole_deg = 0;
for i = 1:length(poly)
    if poly{i}.type == "symmetric"
        poly{i}.p = @(x) (eleSym(x(poly{i}.index),poly{i}.deg));
        poly{i}.grad = @(x) (grad_eleSym(x,poly{i}.index,poly{i}.deg));

    elseif poly{i}.type == "matrix"
        poly{i}.p = @(x) (sum(prod(cat(1,x.^((poly{i}.mat(:,1:end-1)).'),(poly{i}.mat(:,end)).'))));
        poly{i}.deg = sum(poly{i}.mat(1,1:end-1));
        grad_poly = cell(1,n);
        for term = 1:size(poly{i}.mat,1)
            for j = 1:n
                if poly{i}.mat(term,j) ~= 0
                    tmp = poly{i}.mat(term,:); tmp(n+1) = tmp(j)*tmp(n+1); tmp(j) = tmp(j)-1;
                    grad_poly{j}(end+1,:) = tmp;
                end
            end
        end
        poly{i}.grad = @(x) (poly_vec_val(x,grad_poly));
    end
    whole_deg = whole_deg + poly{i}.deg;  %whole_deg is the degree of polynomial
end %whole_deg is the degree of polynomial

mu1 = 1/(12*log(whole_deg)); mu2 = eps/(6*log(whole_deg));

%convert v0 to satisfy minimum eigenvalue of v0 equals to 1/4
eig_v0 = eigH(v0,e,poly);
v0 = e + 0.75*(v0-e)/(1-min(eig_v0));

%%Initialize
itr = 0;
index = 1;
num_grad_cal = 0;
L1 = L/mu1; L2 = L/mu2;
A1 = 0; A2 = 0;
w1 = zeros(n,1); w2 = zeros(n,1); points = zeros(n,1);
w1(:,1) = v0; w2(:,1) = v0; points(:,1) = real(e + (v0-e)/0.75);
v1 = v0; v2 = v0; 
time = zeros(1,1);

%Proj is the orthogonal projection matrix to the hyperplane parallel to {x|Ax=b,<c,x>=<c,v0>}
ker = null(cat(1,c.',A));
Proj = ker*(ker.');

%%main algorithm
while(1)
    %If each of the max iteration or max gradient calculation is attained ,
    %then stop.
    if (itr >= max_itr), fprintf("Reached max iteration \n"), break, end
    if (num_grad_cal > max_num), fprintf("Reached max number of gradient calculation \n"), break, end
    if (time(end) > max_time), fprintf("Reached max time \n"), break, end
    
    %increment the index
    itr = itr + 1;
    
    %update w1_k by accelerated gradient descent
    while(1)
        %set parameters and calculate y
        a = (1+sqrt(1+2*A1*L1))/L1; y = (A1*w1(:,index)+a*v1)/(A1+a);
        
        %calculate eigenvalues of y
        eig_y = eigH(y,e,poly);
        
        %calculate x and increase the number of gradients calculation
        x = y + Proj*real(grad_f_mu(-y,-eig_y,poly,whole_deg,e,mu1,zero_eps))/L1;
        num_grad_cal = num_grad_cal + 1;
        
       %calculate eigenvalues of x
        eig_x = eigH(x,e,poly);
        
        %calculate the gradient at x and increase the number of gradients calculation
        grad_f_hat_mu_x = real(grad_f_mu(-x,-eig_x,poly,whole_deg,e,mu1,zero_eps));
        num_grad_cal = num_grad_cal + 1;
        %If the condition is satisfied, we exit the loop. Otherwise, L1
        %is doubled and the loop continues.
        if dot(grad_f_hat_mu_x,y-x) <= -(norm(Proj*grad_f_hat_mu_x).^2)/L1+zero_eps
            break
        end
        L1 = 2*L1;
    end
    
    %update w1 and reset parameters
    w1(:,index+1) = x; v1 = v1 + a * Proj * grad_f_hat_mu_x; L1 = L1/2; A1 = A1 + a;

    %update w2_k by accelerated gradient descent
    while(1)
        %set parameters and calculate y
        a = (1+sqrt(1+2*A2*L2))/L2; y = (A2*w2(:,index)+a*v2)/(A2+a);
        
        %calculate eigenvalues of y
        eig_y = eigH(y,e,poly);
        
        %calculate x and increase the number of gradients calculation
        x = y + Proj*real(grad_f_mu(-y,-eig_y,poly,whole_deg,e,mu2,zero_eps))/L2;
        num_grad_cal = num_grad_cal + 1;
        
       %calculate eigenvalues of x
        eig_x = eigH(x,e,poly);
        
        %calculate the gradient at x and increase the number of gradients calculation
        grad_f_hat_mu_x = real(grad_f_mu(-x,-eig_x,poly,whole_deg,e,mu2,zero_eps));
        num_grad_cal = num_grad_cal + 1;
        
        %If the condition is satisfied, the exit the loop. Otherwise, L1
        %is doubled and the loop continues.
        if dot(grad_f_hat_mu_x,y-x) <= -(norm(Proj*grad_f_hat_mu_x).^2)/L2 + zero_eps
            break
        end
        L2 = 2*L2;
    end
    
    %update w2 and reset parameters
    w2(:,index+1) = x; v2 = v2 + a * Proj * grad_f_hat_mu_x; L2 = L2/2; A2 = A2 + a;
    
    %calculate points in the original space and objective value
    eig_w2 = eigH(w2(:,index+1),e,poly);
    points(:,index+1) = real(e + (w2(:,index+1)-e)/(1-min(real(eig_w2))));
    time(index+1) = toc;
    
    index = index+1;
    
    %if stopping criteria is satisfied, then stop
    if (f(points(:,index))<threshold), break, end
    
    %calculate eigenvalues of w1_{k+1}
    eig_w1 =  eigH(w1(:,index),e,poly);
    min_eig = min(eig_w1);
    
    %if the minimum eigenvalue of w1_{k+1} is larger than 1/2, conduct the radial update of v_l
    if min_eig >= 0.5-zero_eps
        %fprintf("Radial update occurs at %d iteration.\n",itr)
        w1(:,index+1) = e + 0.75*(w1(:,index)-e)/(1-min_eig);
        w2(:,index+1) = w1(:,index+1);
        points(:,index+1) = real(e + (w2(:,index+1)-e)/0.75);
        time(index+1) = toc;
        index = index + 1;
    end
    %if stopping criteria is satisfied, then stop
    if (f(points(:,index))<threshold), break, end
    
end
end

%% subfunctions
%%%%% subfunctions for the gradient of f_mu %%%%%%
function[grad] = grad_f_mu(x,eigs,poly,whole_deg,e,mu,zero_eps)
%decide m_j
eigs = sort(real(eigs),'descend');
eigs_mult = [eigs(1),1];%Each row corresponds the dissimilar eigenvalue. First column represents the eigenvalue, and second one do its multiplicity. 
l = 1;
for i = 2:length(eigs)
    if (eigs_mult(l,1) - eigs(i)) < zero_eps, eigs_mult(l,2) = eigs_mult(l,2) + 1;
    else eigs_mult(l+1,1) = eigs(i); eigs_mult(l+1,2) = 1; l = l + 1;end
end
mj_exp = exp((eigs_mult(:,1)-eigs_mult(1,1))/mu).*eigs_mult(:,2);
p_mj = zeros(l,1);
grad_p_mj = zeros(length(x),l);
for i = 1:l
    var = x-eigs_mult(i,1)*e;
    p_mj(i,1) = p_deriv(var,eigs_mult(i,2),poly,whole_deg,e);
    grad_p_mj(:,i) = grad_deriv_poly(var,eigs_mult(i,2)-1,poly,whole_deg,e);
end
grad = grad_p_mj*(mj_exp./p_mj)/sum(mj_exp);
end

function[val] = p_deriv(x,r,poly,whole_deg,e)
if int16(r) == 0
    val = 1;
    for i = 1:length(poly)
        val = val * poly{i}.p(x);
    end
else
    tmp = ones(whole_deg,1);
    for j = 1:whole_deg
        for i = 1:length(poly)
            tmp(j) = tmp(j) * poly{i}.p(x+exp(2i*pi*j/whole_deg)*e);
        end
    end
    val = exp(-2i*pi*[1:whole_deg]*r/whole_deg)*tmp*factorial(r)/whole_deg;
end
end