function [model] = AKNS_odemodel(xspan,p,q)
    model=struct('solution',@solution,'Y',@fun_Y,'Z',@fun_Z,'F',@fun_F,'solution_zeros',@solution_zeros);
    function [x,F] = solution(e)
        F0 = [0;1];
        B = [0, -1; 1, 0];
        opts = odeset('Mass', B,'RelTol',1e-8,'RelTol',1e-8 );
        [x,F] = ode45(@(x,F) dFdx(x,F,p,q,e), xspan, F0, opts);
    end
    function [F_end] = fun_F(E)
        F_end = zeros(length(E),2);
        for i = 1:length(E)
            [~,F] = solution(E(i));
            F_end(i,:) = F(end,:);
        end
    end
    function [Y] = fun_Y(E)
        Y = zeros(size(E));
        for i = 1:length(E)
            [~,F] = solution(E(i));
            Y(i) = F(end,1);
        end
    end
    function [Z] = fun_Z(E)
        Z = zeros(size(E));
        for i = 1:length(E)
            [~,F] = solution(E(i));
            Z(i) = F(end,2);
        end
    end
    function [n] = solution_zeros(e)
        F0 = [0;1];
        B = [0, -1; 1, 0];
        opts = odeset('Mass', B,'Events', @Yzeros);
        [~,~,xe] = ode45(@(x,F) dFdx(x,F,p,q,e), xspan, F0, opts);
        n = length(xe);
        function [position,isterminal,direction] = Yzeros(~, F)
            position = F(1);
            isterminal = 0;
            direction = 0;
        end
    end
    function dF = dFdx(x,F,p,q,e)
    dF = [e+q(x), -p(x); -p(x), e-q(x)]*F;
    end
end
