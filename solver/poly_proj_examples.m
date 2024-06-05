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
rng(2024);
n = 50;
deriv_num = 47;
% Generate a point to project
x = randn(n,1);
% format long;
fprintf("Initial point:\n");disp(x);

% Example 1:
% Find the projection using poly_proj
%
% Set up the polynomial.
% In this example this is an elementary symmetric polynomial 
% given in matrix form and corresponds to the derivative relaxation of 
% order deriv_num of the nonnegative orthant.
%
clearvars p
combi = nchoosek(1:n,n-deriv_num);
row = size(combi,1);
p{1}.type = "matrix";
p{1}.mat = zeros(row,n+1);
for i = 1:row
for j = combi(i,:)
    p{1}.mat(i,j) = 1;
end
end
p{1}.mat(:,n+1) = ones(row,1);
e = ones(n,1);

fprintf("Projection:\n");
tic;
[z,status] = poly_proj(x,p,e);
toc;
disp(z);
fprintf("Status:%d\n",status.status);
fprintf("FW Gap:%g\n\n",status.FW_gap);

% Example 2:
% We set the same polynomial but this time we use 
% the elementary symmetric polynomial routines.
% This also illustrates how to pass an optional argument to 
% the solver.
%
clearvars p opts
p{1}.type = "symmetric";
p{1}.deg = n-deriv_num;
p{1}.index = 1:n;
opts.stop.gap_tol = 1e-2; %This is the default FW gap tolerance
%opts.verbose_freq = 1;
e = ones(n,1);
fprintf("Projection using elementary symmetric polynomial routines:\n");
tic;
[z,status] = poly_proj(x,p,e,opts);
toc;
disp(z);
fprintf("Status:%d\n",status.status);
fprintf("FW Gap:%g\n\n",status.FW_gap);

% Example 3:
% Next is a toy example where we compute 
% the projection onto the nonnegative orthant 
% in order to test the oracle interface
%
clearvars p opts
p{1}.type = "oracle";
p{1}.p = @(x) (prod(x));
p{1}.deg = n;
p{1}.grad = @(x) (grad(x));
p{1}.eig = @(x) (min(x));
e = ones(n,1);
fprintf("Projection onto nonnegative orthant:\n");
tic;
[z,status] = poly_proj(x,p,e);
toc;
disp(z);
fprintf("Status:%d\n",status.status);
fprintf("FW Gap:%g\n",status.FW_gap);

function [vec] = grad(x)
    n = length(x);
    vec = zeros(n,1);
    for i = 1:n
        vec(i) = prod(x(setdiff(1:n,i)));
    end
end