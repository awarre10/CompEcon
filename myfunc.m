function [out1,out2,out3] = myfunc(flag,s,x,e,alpha,beta,gamma);
switch flag
    case 'b'; % Bound Function
        out1 = zeros(size(s)); % xl
        out2 = s; % xu
    case 'f'; % Reward Function
        out1 = ((s-x).^(-alpha))/(1-alpha); % f
        out2 = -(s-x).^(-alpha); % fx
        out3 = -alpha*(s-x).^(1-alpha-1); %fxx
    case 'g'; % State Transition Function
        out1 = gamma*x+e.*x.^beta; % g
        out2 = gamma+beta*e.*x.^(beta-1); % gx
        out2 = (beta-1)*beta*e.*x.^(beta-2); % gxx
end