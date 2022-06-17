RK_flag = 4; % Verify order conditions of RK(4, 4)
switch RK_flag
    case 1  % RK(1, 1)
        A = [0; 1]; order = 1;
    case 2  % RK(2, 2), Rlaston's second-order scheme
        A = [0 0; 2/3 0; 1/4 3/4]; order = 2;
    case 21 % RK(2, 2), Heun's second-order scheme
        A = [0 0; 1 0; 1/2 1/2]; order = 2;
    case 3  % RK(3, 3), Heun's third-order scheme
        A = [0 0 0; 1/3 0 0; 0 2/3 0; 1/4 0 3/4]; order = 3;
    case 31 % RK^+(3, 3) in [32]
        alpha = [0 0 0; 1 0 0; 2/3 1/3 0; 59/128+5/128 0 27/64];
        beta = [0 0 0; 2/3 0 0; 0 4/9 0; 15/128*4/3 0 27/64*4/3];
        [A, b] = shuosher2butcher(alpha, beta);
        A = [A; b']; order = 3;
    case 32 % RK(3, 3) in [54]
        A = [0 0 0; 1 0 0; 1/4 1/4 0; 1/6 1/6 2/3]; order = 3;
    case 4  % RK(4, 4)
        A = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0; 1/6 1/3 1/3 1/6]; order = 4;
    case 41 % RK(5, 4) in [39]
        A = [0 0 0 0 0;2/5 0 0 0 0; 1/10 1/2 0 0 0; 1/16 1/16 3/8 0 0; ...
            1/10 1/10 4/15 8/15 0; 5/32 25/96 25/96 1/6 5/32]; order = 4;
    case 42 % RK^+(5,4) in [32]
         alpha = [0 0 0 0 0;
             0.387392167970373+0.612607832029627 0 0 0 0;
             0.568702484115635 0.431297515884365 0 0 0;
             0.589791736452092 0 0.410208263547908 0 0;
             0.213474206786188 0 0 0.786525793213812 0;
             0.270147144537063+0.029337521506634 0.239419175840559 0 ...
             0.227000995504038 0.234095162611706];
         r = 1.346586417284006;
         beta = [0 0 0 0 0; 0.612607832029627  0 0 0 0;
             0 0.431297515884365 0 0 0;
             0 0 0.410208263547908 0 0;
             0 0 0 0.786525793213812 0;
             0.029337521506634 0.239419175840559 0 ...
             0.227000995504038 0.234095162611706]/r;
         [A, b] = shuosher2butcher(alpha, beta);
         % The file shuosher2butcher.m can be downloaded via https://www.cfm.brown.edu/people/sg/SSPpage/sspsite/matlab_scripts.html
         A = [A; b']; order = 4;
    case 43 % RK(5, 4) in [57]
        alpha = [0 0 0 0 0;
            1, 0 0 0 0;
            0.444370493651235 0.555629506348765 0 0 0;
            0.620101851488403 0 0.379898148511597 0 0;
            0.178079954393132 0 0 0.821920045606868 0;
            0 0 0.517231671970585 0.096059710526147 0.386708617503268];
        beta = [0 0 0 0 0;
            0.39175222657189  0 0 0 0;
            0 0.368410593050371 0 0 0;
            0 0 0.251891774271694 0 0;
            0 0 0 0.54497475022852  0;
            0 0 0 0.063692468666290 0.226007483236906];
        [A, b] = shuosher2butcher(alpha, beta);
        A = [A; b']; order = 4;
    case 44 % eSSPRK(10, 4) in [32]
        A = zeros(10, 10);
        for i = 1:10
            for j = 1:10
                if i < 6
                    if j < i
                        A(i,j) = 1/6;
                    end
                elseif i >= 6
                    if j <= 5
                        A(i,j) = 1/15;
                    end
                end
                if j>= 6 && i <=10 && i>= 6 && j < i
                    A(i,j) = 1/6;
                end
            end
        end
        b = 1/10*ones(10, 1);
        A = [A; b']; order = 4;
    otherwise
        stage = 4; order = 4;
        A = sym('A', [stage+1, stage]);
        for i = 1:stage
            for j = i:stage
                A(i,j) = 0;
            end
        end
end
syms x;
stage = size(A, 2);
A_hat = sym(zeros(size(A)));
c_hat = sym(zeros(stage+1,1));
psi   = sym(zeros(stage+1,1));

% Calculate coefficients of pRK method
psi(1) = 1;
for i = 2:stage+1
    psi(i) = 1; c_hat(i) = 0;
    for j = 1:i-1
        psi(i) = psi(i) + x * A(i,j) * psi(j);
    end
    for j = 1:i-1
        A_hat(i,j) = A(i,j) * psi(j) / psi(i);
        for k = j+1:i-1
            A_hat(i,j) = A_hat(i,j) + x * A(i,k) * psi(k) / psi(i) * A_hat(k,j);
        end
        c_hat(i) = c_hat(i) + A_hat(i,j);
    end
end

% Verification of order conditions
b_hat = A_hat(stage+1,:);
A_hat = A_hat(1:stage,:);
c_hat = c_hat(1:stage);
c     = sum(A(1:stage,:), 2);
if order >= 1
    fprintf('Leading Error x^%d:\n', order);
    taylor(sum(b_hat) - 1, x, 0, 'order', order+2)
end
if order >= 2
    fprintf('Leading Error x^%d:\n', order-1);
    for i = 0:1
        taylor(b_hat*(c_hat.^i.*c.^(1-i)) - 1/2, x, 0, 'order', order+1)
    end
end
if order >= 3
    fprintf('Leading Error x^%d:\n', order-2);
    for i = 0:2
        taylor(b_hat*(c_hat.^i.*c.^(2-i)) - 1/3, x, 0, 'order', order)
    end
    for i = 0:1
        taylor(b_hat*A_hat*(c_hat.^i.*c.^(1-i)) - 1/6, x, 0, 'order', order)
    end
end
if order >= 4
    fprintf('Leading Error x^%d:\n', order-3);
    for i = 0:3
        taylor(b_hat*(c_hat.^i.*c.^(3-i)) - 1/4, x, 0, 'order', order-1)
    end
    for i = 0:1
        for j = 0:1
            taylor(b_hat*((c_hat.^i.*c.^(1-i)).*(A_hat*(c_hat.^j.*c.^(1-j)))) - 1/8, x, 0, 'order', order-1)
        end
    end
    for i = 0:2
        taylor(b_hat*A_hat*(c_hat.^i.*c.^(2-i)) - 1/12, x, 0, 'order', order-1)
    end
    for i = 0:1
        taylor(b_hat*A_hat^2*(c_hat.^i.*c.^(1-i)) - 1/24, x, 0, 'order', order-1)
    end
end
