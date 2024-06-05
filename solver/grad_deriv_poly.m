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
function[grad] = grad_deriv_poly(x,mult,poly,degree,e)
if int16(mult) <= 0
    grad =  grad_prod_poly(x,poly);
else
    scaler = 1e5/max(abs(x));
    x = x*scaler;
    tmp_mat = zeros(length(x),degree);
    for j = 1:degree
        tmp_mat(:,j) = grad_prod_poly(x+exp(2i*pi*j/degree)*e,poly);
    end
    grad = tmp_mat*exp(-2i*pi*[1:degree].'*mult/degree);
    grad = grad*factorial(mult)/(degree*(scaler.^(degree-mult-1)));
end
end

function[grad] = grad_prod_poly(x,poly)
%grad_poly{i} is the polynomial vector representing the gradient of poly{i}
%grad is grad(p_1*...*p_l) = grad(p_1)*p_2*...*p_l + p_1*grad(p_2)*...*p_l + ...  

poly_vals = zeros(1,length(poly));
for i = [1:length(poly)]
    poly_vals(i) = poly{i}.p(x);
end
pr = prod(poly_vals);

grad = zeros(size(x));
for i = 1:length(poly)
    if isfinite(pr/poly_vals(i))
        grad = grad + poly{i}.grad(x)*pr/poly_vals(i);
    else
        grad = grad + poly{i}.grad(x)*prod(poly_vals([1:i-1,i+1:end]));
    end
end
end