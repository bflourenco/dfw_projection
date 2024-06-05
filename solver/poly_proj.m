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
function[x,status] = poly_proj(x0,p,e,opts)
%Computes a low-to-medium accuracy orthogonal projection of x0 onto the hyperbolicity 
%cone given by p and e using the dual Frank-Wolfe method in FW_HP.m.
%By default, x0 is scaled to have 2-norm less or equal than 1 before being passed 
%to the solver and, then, the output is unscaled back. This can 
%be disabled by setting opts.scale = False.
% 
%See poly_proj_examples.m for examples.
%
%Input:
%x0 - [m by 1 vector]
%p - set of polynomials following the format described in FW_HP.m
%e - [m by 1 vector]
%opts - options to be passed to the solver, see opts argument in FW_HP.m.
%Addionally, it is possible to set opts.scale to false in order to disable 
%input scaling.
%
%
%Output:
%x - [m by 1 vector] 
%status: it is the same as the output argument "status" of FW_HP.m:
%       status - 1: Obtained a feasible solution with FW gap <= stop.gap_tol
%                2: Obtained a feasible solution but FW gap stop criterion
%                was not satisfied. In this case, x is the feasible solution 
%                found with the smallest FW gap.
%                3: No primal feasible solution found. In this case, 
%                x is the solution found with the smallest feasibility
%                violation.
%       time: total running time.
%       itr: total number of iterations.
%       FW_gap: Best FW_gap found. It is Inf if no dual feasible solution
%       was found
%

n = size(x0,1);
if nargin < 4
    opts = {};
end
if isfield(opts,'step') == false, opts.step = {}; end
if isfield(opts.step,'rule') == false, opts.step.rule = "Lipschitz"; end
if isfield(opts.step,'L') == false, opts.step.L = 1; end
if isfield(opts,'stop') == false, opts.stop = {}; end
if isfield(opts.stop,'max_time') == false, opts.stop.max_time = 10; end
if isfield(opts,'y0') == false, opts.y0 = zeros(n,1); end   
if isfield(opts,'scale') == false, opts.scale = true; end   

scale_factor = 1;
if opts.scale == true && norm(x0) > 1  
    scale_factor = norm(x0);
end
x1 = x0/scale_factor;


grad = @(x) (x+x1);
c = norm(e)*sqrt(norm(e-x1).^2+norm(e).^2);
if dot(e,opts.y0) > c, c = dot(e,opts.y0); end

[x,~,status] = FW_HP(grad,eye(n),zeros(n,1),p,e,c,opts);
x = x*scale_factor;

end

