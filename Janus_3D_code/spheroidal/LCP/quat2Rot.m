function rotConfig = quat2Rot(quatConfig, a, b, c);
    %Converts a vector of translation points and quaternions to a struct of centers, rotations matrices, and semiaxes.
    
    C = cell(length(quatConfig)/7, 1);
    R = cell(length(quatConfig)/7, 1);
    for i = 1:length(quatConfig)/7
        C{i} = quatConfig(7*(i-1)+1:7*(i-1)+3);
        R{i} = quat2rotm(quatConfig(7*(i-1)+4:7*(i-1)+7).');
    end
    rotConfig = struct('a', a, 'b', b, 'c', c,'C', C, 'R', R);
end