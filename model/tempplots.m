h = 15;
W = conv2(R, ones(h,1)/h, 'same');
W = W(1:h:end,:);
C = conv2(B, ones(h,1)/h, 'same');
C = C(1:h:end,:);

U = T(1:h:end,:);

for i = 1:3
    subplot(3,1,i)
    plot(U, W(:, (1:3)+(i-1)*3) - repmat(mean(C,2),1,3))
    ylim([-5,10]);
    xlim([1, 23]);
end
%subplot(1,4,4)
figure
plot(U, C, 'k'), xlim([1, 23]);
