function quatConfig = rot2Quat(rotConfig)
    %Converts a struct of centers, rotations matrices, and semiaxes to a vector of translation points and quaternions.
    quatConfig = zeros(7*length(rotConfig),1);
    for i = 1:length(rotConfig)
        quatConfig(7*(i-1)+1:7*(i-1)+3) = rotConfig(i).C;
        quatConfig(7*(i-1)+4:7*(i-1)+7) = rotm2quat(rotConfig(i).R).';
    end

end