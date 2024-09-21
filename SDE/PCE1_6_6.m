clear all

p = @(lamb) [lamb/2 lamb/2 1-lamb];
P = [.5 .5 0;.5 .5 0;0 0 1];

n = 10000;

for i = 1:n
    Pe(i,:) = p(rand())*P;
end

fprintf('Expected:\n');
disp(p(0.5));
fprintf('Obtained after averaged %d random vector:\n',n);
disp(sum(Pe)/n);
