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
function[a] = stepsize(step,x,grad,descent,itr)
% This function calculates step size. "stepsize_rule.rule" specifies which rule are used. Available step size rules and
% necessary additional arguments are below. Default rule is "deminishing".
% "diminishing": No additional arguments are needed. Step size at k-th iteration is 2/k+2.
% "constant": Constant number used for step size is in "step.stepsize". This variable must be in [0,1]. Might not work if problem 
% data is not strongly convex.
% "optimal": step.stepfunc is the function(x,grad,descent) which calculates optimal
%            step size at x along descent. This function must return the stepsize within [0,1].
% "Armijo": step.f is the objective function f(x).
% "Lipschitz": step.L is a Lipschitz constant of the gradient of the objective function.

switch step.rule
    case "constant"
        a = step.stepsize;
    case "optimal"
        a = step.stepfunc(x,grad,descent);
    case "Armijo"
        a = 1; xi = 1e-3; tau = 0.5; %xi and tau are the hyperparameters
        fx = step.f(x);
        tmp = dot(grad,descent);
        while(step.f(x+a*descent) > fx + xi*a*tmp)
            a = a * tau;
        end
    case "Lipschitz"
        a = -dot(grad,descent)/(step.L*(norm(descent).^2));
    otherwise
        a = 2/(itr+2);
end
    