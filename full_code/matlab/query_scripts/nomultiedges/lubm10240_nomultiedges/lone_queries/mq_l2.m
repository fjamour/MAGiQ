function R = mq_l2(I, J, V, N, GPU, varargin)
t_layers = tic;
[hI3, hJ3] = make_predicate_graph(I, J, V, 3);
[hI6, hJ6] = make_predicate_graph(I, J, V, 6);
fprintf('Time to make layers: %f\n', toc(t_layers));

use_gpu = ismember('use_gpu', varargin);
if use_gpu
    reset(GPU)
    t_gpu_copy = tic;
    I3 = gpuArray(hI3); J3 = gpuArray(hJ3);
    I6 = gpuArray(hI6); J6 = gpuArray(hJ6);
    fprintf('Time to copy layers to GPU: %f\n', toc(t_gpu_copy));
else
    t_copy = tic;
    I3 = hI3; J3 = hJ3;
    I6 = hI6; J6 = hJ6;
    fprintf('Time to make local layers copies: %f\n', toc(t_copy));
end

t_all = tic;
t_mxv = tic;
fv0 = 80;
fv1 = mxv(I6, J6, fv0);
fv2 = mxv(I3, J3, fv1, 'trans');
fv1 = vandv(fv1, mxv(I3, J3, fv2));
fv0 = vandv(fv0, mxv(I6, J6, fv1, 'trans'));
fprintf('Time to do all mxv and vandv: %f\n', toc(t_mxv));

t_bm = tic;
sB01 = make_bm(I6, J6, fv0, fv1, 'trans');
sB12 = make_bm(I3, J3, fv1, fv2);
fprintf('Time to do make binding mats: %f\n', toc(t_bm));

t_R = tic;
R = sB01;
R = R_join(R, sB12, 2);
fprintf('Time to do joins to get R: %f\n', toc(t_R));
fprintf('Total query time: %f\n', toc(t_all));
fprintf('Number of results: %d\n', size(R, 1));

end
