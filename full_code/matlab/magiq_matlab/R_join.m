% expands R by joining R (column R(:, jv)) with BM (column BM(:, 1))
function R_out = R_join(R, BM, jv)
    R = sortrows(R, [jv]);
    [ia, ib] = merge_join( gather(R(:, jv)), gather(BM(:, 1)) );
    new_col = BM(:, 2);
    R_out = [R(ia, :), new_col(ib)];
end
