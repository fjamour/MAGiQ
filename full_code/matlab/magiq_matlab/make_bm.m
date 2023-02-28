% returns sorted (on first dimension) binding matrix for vr and vc
function BM = make_bm(I, J, vr, vc, varargin)
    trans = ismember('trans', varargin);
    if trans % diag(vr) x m' x diag(vr)
        idx = ismember(J, vr) & ismember(I, vc);
        BM  = sortrows([J(idx), I(idx)]);
    else     % diag(vr) x m  x diag(vr)
        idx = ismember(I, vr) & ismember(J, vc);
        BM  = sortrows([I(idx), J(idx)]);
    end
end