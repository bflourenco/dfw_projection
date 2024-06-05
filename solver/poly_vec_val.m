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
function[vec] = poly_vec_val(x,poly_vec)
%vec is [poly_vec{1} at x; ... ;poly_vec{end} at x]
vec = zeros(length(poly_vec),1);
for i = 1:length(poly_vec)
    if isempty(poly_vec{i})
        vec(i) = 0;
    else
        vec(i) = sum(prod(cat(1,x.^((poly_vec{i}(:,1:end-1)).'),(poly_vec{i}(:,end)).')));
    end
end