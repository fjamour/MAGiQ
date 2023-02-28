function v_out = vandv(v1, v2)
    v_out = v1(ismember(v1, v2));
end