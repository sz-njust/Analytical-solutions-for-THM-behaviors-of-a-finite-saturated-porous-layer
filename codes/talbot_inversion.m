function ilt = talbot_inversion(f_s, t, M)

    if size(t, 1) == 1
        t = t';
    elseif size(t, 2) > 1
        error('Input times, t, must be a vector.');
    end

    if nargin < 3
        M = 64;
    end
        
    k = 1:(M-1); 

    delta = zeros(1, M);
    delta(1) = 2*M/5;
    delta(2:end) = 2*pi/5 * k .* (cot(pi/M*k)+1i);

    gamma = zeros(1, M);
    gamma(1) = 0.5*exp(delta(1));
    gamma(2:end) = (1 + 1i*pi/M*k.*(1+cot(pi/M*k).^2)-1i*cot(pi/M*k))...
                   .* exp(delta(2:end));
    
    [delta_mesh, t_mesh] = meshgrid(delta, t);
    gamma_mesh = meshgrid(gamma, t);
    
    ilt = 0.4./t .* sum(real(   gamma_mesh ...
                             .* arrayfun(f_s, delta_mesh./t_mesh)), 2);

end