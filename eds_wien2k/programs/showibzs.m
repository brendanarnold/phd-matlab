function showIBZs(IBZklist,rlv,symmops);

figure;

cols=([1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 0.5 0.5 0; 0.5 0 0.5; 0 0.5 0.5; 1 0.5 0; 0.5 1 0; 0.5 0 1]);
hold on;
for n=1:size(symmops,3)
    symmIBZ{n}=IBZklist*symmops(:,:,n)';
    cartIBZ{n}=[symmIBZ{n}(:,1)*rlv(1,1)+symmIBZ{n}(:,2)*rlv(2,1) ...
        symmIBZ{n}(:,1)*rlv(1,2)+symmIBZ{n}(:,2)*rlv(2,2)]/28;
    plot(cartIBZ{n}(:,1),cartIBZ{n}(:,2),'o','MarkerEdgeColor',cols(n,:));
end



