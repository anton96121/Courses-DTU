function[tau] = R(lambda,h,b,A)
    tau = (1+lambda*h*(b)'*inv(eye(length(b))-lambda*h*A)*ones(length(b),1));
end