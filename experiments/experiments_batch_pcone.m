tStart = tic;

n = [100,300,500,1000];
p = [1.1, 1.3, 3,5];
num_points = 30;
ratio_list = [1e-1,5e-2,1e-2,5e-3,1e-3];
skip_DDS = true;

for i = 1:length(n)
    for j = 1:length(p)
        ex_proj_pcone_table(n(i),p(j),num_points,1e-8,ratio_list,skip_DDS,true);
        ex_proj_pcone_table(n(i),p(j),num_points,1e-3,ratio_list,skip_DDS,true);
    end
end
% ex_proj_pcone_table(50,1.1,30,1e-8,[1e-1,5e-2,1e-2,5e-3,1e-3],true,true);
% ex_proj_pcone_table(50,1.1,30,1e-3,[1e-1,5e-2,1e-2,5e-3,1e-3],true,true);
% ex_proj_pcone_table(100,1.1,30,1e-8,[1e-1,5e-2,1e-2,5e-3,1e-3],true,true);
% ex_proj_pcone_table(100,1.1,30,1e-3,[1e-1,5e-2,1e-2,5e-3,1e-3],true,true);
% ex_proj_pcone_table(250,1.1,30,1e-8,[1e-1,5e-2,1e-2,5e-3,1e-3],true,true);
% ex_proj_pcone_table(250,1.1,30,1e-3,[1e-1,5e-2,1e-2,5e-3,1e-3],true,true);
% ex_proj_pcone_table(500,1.1,30,1e-8,[1e-1,5e-2,1e-2,5e-3,1e-3],true,true);
% ex_proj_pcone_table(500,1.1,30,1e-3,[1e-1,5e-2,1e-2,5e-3,1e-3],true,true);
% ex_proj_pcone_table(1000,1.1,30,1e-8,[1e-1,5e-2,1e-2,5e-3,1e-3],true,true);
% ex_proj_pcone_table(1000,1.1,30,1e-3,[1e-1,5e-2,1e-2,5e-3,1e-3],true,true);

toc(tStart);