function [V,Vn,J] = JCF(A)

% clc
% close all
% clear all
%
% A_ideal = [-1 2 0 0 1 0;
%     -2 -1 0 0 0 1;
%     0 0 -1 1 0 0;
%     0 0 0 -1 0 0;
%     0 0 0 0 -1 2;
%     0 0 0 0 -2 -1]
% T = randn(length(A_ideal));
% A = T*A_ideal/T;

n = length(A);
[Vm,Jm] = jordan(A);
temp_eig = diag(Jm);

%%
k = 1;
Lia = 1:1:n;
while k <= n
    index{k} = [];
    for j = 1:length(Lia)
        [V, li] = AddVectors(Vm(:,min(Lia)),Vm(:,Lia(j)));
        if li == 0
            index{k} = [index{k} Lia(j)];
        end
    end
    ind = [];
    for i = 1:k
        ind = [ind index{i}];
    end
    ind = sort(ind);
    Lia = find(not(ismember([1:1:n],ind)));
    if isempty(Lia)
        break
    end
    k = k+1;
end

%%
J = [];
V = [];
for i = 1:length(index)
    ind = index{i};
    q = length(index{i});
    if isreal(temp_eig(index{i}))
        M = diag(temp_eig(index{i}))+diag(ones(q-1,1),1);
        for j = 1:q
            if j == 1
                V = [V Vm(:,ind(j))];
            else
                temp = pinv(A-temp_eig(ind(1))*eye(n))*V(:,end);
                V = [V temp];
            end
        end
    else
        M = [];
        for j = 1:2:q
            alpha = real(temp_eig(ind(j)));
            beta  = imag(temp_eig(ind(j)));
            M = blkdiag(M,[alpha beta; -beta alpha]);
        end
        M = M +kron(diag(ones(q/2-1,1),1),eye(q/2));
        Vtemp = [];
        for j = 1:q/2
            if j == 1
                Vtemp = [Vtemp Vm(:,ind(j))];
            else
                temp = pinv(A-temp_eig(ind(1))*eye(n))*Vtemp(:,end);
                Vtemp = [Vtemp temp];
            end
        end
        V = [V MakeReal(Vtemp)];
    end
    J = blkdiag(J,M);
end
%J
%V
%test = V\A*V

Vtemp = MakeReal(Vm);
if rank(Vtemp)== n && rcond(Vtemp) > 1e-4
    V = Vtemp;
    J = V\A*V;
end

Vn = Normalisation(V);

 end

%% SUPPORT FUNCTIONS

function [res, li] = AddVectors(V, v_candidate)

if norm(imag(v_candidate)) < 1e-4
    V_test = [V v_candidate];
else
    V_test = [V real(v_candidate) imag(v_candidate)];
end

svd_test = svd(V_test);
rcond_test = min(svd_test)/max(svd_test);
if rcond_test < 1e-4
    V_test = V;
end

size_test = size(null(V_test.'));
if isempty(V)
    size_M = length(v_candidate)*[1 1];
else
    size_M = size(null(V.'));
end
if size_test(2) < size_M(2)
    res = V_test;
    li = 1;
else
    res = V;
    li = 0;
end

end
%
% function res = GeneralisedEigVector(M,v)
% N = null(M);
% [n m] = size(N);
% T = [N null(N.')];
% barM = T\M*T;
% barv = T\v;
% v_candidate = barM(1:m,m+1:end)*null(barM(m+1:end,m+1:end));
% if isempty(v_candidate)
%     res = [];
% else
%     test =  (v_candidate./norm(v_candidate)).'*(barv(1:m)/norm(barv(1:m)));
%     if abs(test) < 1-1e-4
%         res = [];
%     else
%         res = T\[zeros(m,1); v_candidate]*(norm(barv(1:m))/norm(v_candidate));
%     end
% end
%
% end

function V_real = MakeReal(V)
V_real = V*0;
n = size(V);
k = 1;
for i = 1:n(2)
    if sum(abs(imag(V(:,i)))) < 1e-4
        V_real(:,k) = real(V(:,i));
        k = k+1;
    else
        V_real(:,k) = real(V(:,i));
        V_real(:,k+1) = imag(V(:,i));
        k = k+2;
    end
end
size_V_real = size(V_real);
if size_V_real(2) > n(1)
    V_real = V_real(:,1:n(1));
end
end

function Vn = Normalisation(V)
Vn = V*0;
for kk = 1:length(Vn)
    nV = norm(V(kk,:),1);
    Vn(kk,:) = abs(V(kk,:)./nV)*100;
end
end