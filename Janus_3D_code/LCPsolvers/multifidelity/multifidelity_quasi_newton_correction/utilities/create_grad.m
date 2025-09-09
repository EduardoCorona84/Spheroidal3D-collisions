function grad = create_grad(fg)
    %this function undoes create_fg and only gives back the gradient (primarily useful for the inner solver)
    grad = @temp_grad;
    function g = temp_grad(x_k)
        [~, g, ~] = fg(x_k, [], [], []);

    end
end