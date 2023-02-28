function R = mq_l1(I, J, V, N, GPU, varargin)
t_layers = tic;
[hI2, hJ2] = make_predicate_graph(I, J, V, 2);
[hI6, hJ6] = make_predicate_graph(I, J, V, 6);
[hI10, hJ10] = make_predicate_graph(I, J, V, 10);
[hI11, hJ11] = make_predicate_graph(I, J, V, 11);
fprintf('Time to make layers: %f\n', toc(t_layers));

use_gpu = ismember('use_gpu', varargin);
if use_gpu
    reset(GPU)
    t_gpu_copy = tic;
    I2 = gpuArray(hI2); J2 = gpuArray(hJ2);
    I6 = gpuArray(hI6); J6 = gpuArray(hJ6);
    I10 = gpuArray(hI10); J10 = gpuArray(hJ10);
    I11 = gpuArray(hI11); J11 = gpuArray(hJ11);
    fprintf('Time to copy layers to GPU: %f\n', toc(t_gpu_copy));
else
    t_copy = tic;
    I2 = hI2; J2 = hJ2;
    I6 = hI6; J6 = hJ6;
    I10 = hI10; J10 = hJ10;
    I11 = hI11; J11 = hJ11;
    fprintf('Time to make local layers copies: %f\n', toc(t_copy));
end

t_all = tic;
t_mxv = tic;
fv0 = 22639;
fv5 = mxv(I6, J6, fv0);
fv3 = mxv(I10, J10, fv5);
fv1 = 25;
fv3 = vandv(fv3, mxv(I6, J6, fv1));
fv4 = mxv(I2, J2, fv3, 'trans');
fv2 = 8622223;
fv4 = vandv(fv4, mxv(I6, J6, fv2));
fv6 = mxv(I11, J11, fv4);
fv5 = fv6;
fv4 = vandv(fv4, mxv(I11, J11, fv6, 'trans'));
fv3 = vandv(fv3, mxv(I2, J2, fv4));
fv5 = vandv(fv5, mxv(I10, J10, fv3, 'trans'));
fv0 = vandv(fv0, mxv(I6, J6, fv5, 'trans'));
fprintf('Time to do all mxv and vandv: %f\n', toc(t_mxv));

t_bm = tic;
sB05 = make_bm(I6, J6, fv0, fv5, 'trans');
sB53 = make_bm(I10, J10, fv5, fv3, 'trans');
sB34 = make_bm(I2, J2, fv3, fv4);
sB46 = make_bm(I11, J11, fv4, fv6, 'trans');
sB42 = make_bm(I6, J6, fv4, fv2);
sB31 = make_bm(I6, J6, fv3, fv1);
fprintf('Time to do make binding mats: %f\n', toc(t_bm));

t_R = tic;
R = sB05;
R = R_join(R, sB53, 2);
R = R_join(R, sB34, 3);
R = R_filter(R, sB46, 4, 2);
R = R_join(R, sB42, 4);
R = R_join(R, sB31, 3);
fprintf('Time to do joins to get R: %f\n', toc(t_R));
fprintf('Total query time: %f\n', toc(t_all));
fprintf('Number of results: %d\n', size(R, 1));

end
