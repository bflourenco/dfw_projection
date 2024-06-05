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
function [vec] = grad_eleSym(x,indexes,k)
%grad_eleSym Return the gradient of the k-th elementary symmetric polynomial of x(indexes).
%This function uses the fact that 
% d eleSym(x(indexes),k) /d x_i = 0 if i is not in indexes, and
% d eleSym((x_1,...,x_n),k) /d x_i = eleSym((x_1,...,x_i-1,x_i+1,...,x_n),k-1) otherwise.
n = length(x);
vec = zeros(n,1);
for i = 1:n
    if ismember(i,indexes)
        vec(i) = eleSym(x(setdiff(indexes,i)),k-1);
    end
end
end 