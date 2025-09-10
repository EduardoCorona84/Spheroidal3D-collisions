function A = create_A(fg)
    %this function undoes create_fg and only gives back the matvec

    A = @temp_A;
    function Ax = temp_A(x_k)
        [~, ~, Ax] = fg(x_k, [], [], []);

    end
end