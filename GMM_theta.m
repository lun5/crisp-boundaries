AIC = zeros(1,3);
BIC = zeros(1,3);
obj = cell(1,3);
options = statset('Display','final','MaxIter',200);
for k = 1:3
    obj{k} = gmdistribution.fit(theta',k,'Options',options);
    AIC(k)= obj{k}.AIC;
    BIC(k) = obj{k}.BIC;
end

[minAIC,numComponents_AIC] = min(AIC);
[minBIC,numComponents_BIC] = min(BIC);
numComponents = max(numComponents_AIC,numComponents_BIC);
numComponents

model = obj{numComponents};
mus = model.mu;
sigmas = model.Sigma;
disp('Means are'); mus(:)
disp('Standard devs are'); sigmas(:)

figure;
col_vec = {'r','g','b'};
binranges = -pi:0.1:pi;
[bincounts] = histc(theta,binranges);
figure;bar(binranges,bincounts/sum(bincounts),'histc');
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0.8 .8 .8],'EdgeColor','w');
xlabel('Theta'); 
      
hold on;
y = cell(1,3);
x = -pi:0.1:pi;
for i = 1:3
    y{i} = normpdf(x,mus(i),sigmas(i));
    plot(x,y{i}./length(x),col_vec{i},'LineWidth',3); hold on;
end
hold off;

p = kde(theta,0.05,[],'e');
pd = evaluate(p,X(:)',tol);
figure; 
bar(binranges,bincounts/sum(bincounts),'histc');
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0.8 .8 .8],'EdgeColor','w');
xlabel('Theta');
hold on
plot(X(:),pd./10,'LineWidth',3);
hold off
