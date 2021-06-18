function [fun_zeros] = find_zeros(fun,xspan)
X = xspan(1):xspan(2)+1;
fun_sign = sign(fun(X));
leftshift_fun_sign = [fun_sign(2:end),0];
index_change_fun_sign = find(fun_sign.*leftshift_fun_sign<0);
fun_zeros = zeros(length(index_change_fun_sign),1);
for k = 1:length(index_change_fun_sign)
    index = index_change_fun_sign(k);
    fun_zeros(k) = fzero(fun,[X(index),X(index+1)]);
end
end

