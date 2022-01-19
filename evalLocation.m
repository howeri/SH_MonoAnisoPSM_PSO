function cost = evalLocation2(location)
    location = location';
    global H
    global f
    global direction
    global res_Re
    global res_Im
    warning('off','all')
    
    %% Original Evaluation
%     tic;
%     lead = 2*direction*location/physconst('LightSpeed');
%     a = exp(1i*2*pi*lead*f);
%     e = [reshape(res_Re,[14641*size(f, 2),1]); reshape(res_Im,[14641*size(f, 2),1]) ];
%     S = repmat(H, size(f, 2)*2, 1).*[real(reshape(a, [14641*size(f, 2),1])) ; imag(reshape(a, [14641*size(f, 2),1]))];
%     x_hat = S\e;
%     cost = norm(S*x_hat-e);
%     toc;
    
    %% Efficient Evaluation
    lead = 2*direction*location/physconst('LightSpeed');
    w = exp(1i*2*pi*lead*f);
    U = size(f,2)*(H'*H);
    T = U\H';
    alpha_n = zeros(size(U,1),1);
    for i=1:size(f, 2)
        alpha_n = alpha_n + T*(real(w(:,i)).*res_Re(:,i));
        alpha_n = alpha_n + T*(imag(w(:,i)).*res_Im(:,i));
    end
     
    cost=0;
    for i=1:size(f, 2)
        cost = cost + sum((res_Re(:,i) - (real(w(:,i)).*H)*alpha_n).^2);
        cost = cost + sum((res_Im(:,i) - (imag(w(:,i)).*H)*alpha_n).^2);
    end
    cost = sqrt(cost); 
end

