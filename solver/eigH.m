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
function[vec] = eigH(x,e,poly)
%This function is a subfunction to calculate eigenvalues of x for (poly,e).
vec = [];
for i = 1:length(poly)
    if isfield(poly{i},"eig")
        vec = [vec;real(poly{i}.eig(x))];
    else
        if poly{i}.type == "symmetric"
            vec = [vec;real(eig_eleSym(x,e,poly{i}))];
        else
            vec = [vec;real(eig_general(x,e,poly{i}))];
        end
    end
end
end

function[eig] = eig_general(x,e,poly)
if poly.deg == 1;
    eig = [poly.p(x)/poly.p(e)];
else
    %scale x
    scaler = norm(x,Inf);
    x = x / scaler;

    %calculate the coefficients of p(x-te) by inverse Fourier transform
    values = zeros(poly.deg,1);
    f0 = poly.p(x);
    for j = 1:poly.deg
        values(j) = poly.p(x-exp(2i*(j-1)*pi/poly.deg)*e)-f0;
    end
    coefficients = [real(ifft(values));f0];

    %the roots of p(x-te) are the eigenvalues of x
    eig = roots(coefficients)*scaler;
end
end

function[eig] = eig_eleSym(x,e,poly)
if poly.deg == 1
    eig = [poly.p(x)/poly.p(e)];
else
    x = x(poly.index);
    %scale x
    scaler = norm(x,1);
    x = x / scaler;

    %calculate the coefficients of p(x+te) by convolution
    n = length(x);
    queue = cell(1,n);
    for i = 1:n
        queue{i} = [x(i) 1];
    end
    while size(queue,2) > 1
        p = conv(queue{1},queue{2});
        deg_p = size(p,2) - 1;
        if deg_p > poly.deg
            queue = {queue{3:end},p((deg_p-poly.deg+1):end)};
        else
            queue = {queue{3:end},p};
        end
    end
    coefficients = flip(queue{1});
    c = 1; tmp = n-poly.deg;
    for i = 1:poly.deg
        c = c*(tmp+i)/i;
        coefficients(end-i) = c*coefficients(end-i);
    end
    %minus the roots of p(x+te) are the eigenvalues of x
    eig = -roots(coefficients)*scaler;
end
end