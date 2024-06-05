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
function [value] = eleSym(x,k)
%ELESYM Return the value of the k-th elementary symmetric polynomial at x
%
% The idea is that the value of the k-th elementary symmetric polynomial at
% x corresponds to coefficient of y^k in \prod_{i=1}^n (1+x_i*y). Then,
% we evaluate \prod_{i=1}^n (1+x_i*y) through a
% simple divide-and-conquer algorithm.
% Polynomial multiplication is done with the "conv"
% function. An approach using FFTs would have a better complexity, but
% in our experiments, unless the degree is very large,  "conv" is faster.
    if isrow(x)
        n = size(x,2);
    else
        n = size(x,1);
    end
    queue = cell(1,n);
    for i = 1:n
        queue{i} = [x(i) 1];
    end
    while size(queue,2) > 1
        p = conv(queue{1},queue{2});
        deg_p = size(p,2) - 1;
        if deg_p > k
            queue = {queue{3:end},p((deg_p-k+1):end)};
        else
            queue = {queue{3:end},p};
        end
    end
    value = queue{1}(1);
end