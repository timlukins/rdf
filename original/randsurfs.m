function randsurfs 

s = 0.05;
[x,y] = meshgrid(-3:s:3,-3:s:3);
z = peaks(x,y);

[I,II]=fundforms(x,y,z,s);

%zi = peaks(x,y);
zi = peaks(x,y)*1.1;
%zi = z';

[Ii,IIi]=fundforms(x,y,zi,s);

for u=1:size(I,1)
  for v=1:size(I,1)

    A(:,:) = Ii(u,v,:,:)  - I(u,v,:,:);
    B(:,:) = IIi(u,v,:,:) - II(u,v,:,:);
    
    Da(u,v) = rank(A);
    Db(u,v) = rank(B);

  end;
end;

figure; subplot(2,1,1); imshow(Da,[1 3]); hold on; colormap jet(3); colorbar; title('A rank'); 
subplot(2,1,2); imshow(Db,[1 3]); hold on; colormap jet(3); colorbar; title('B rank');
%figure; hold on; surf(z,'CData',zeros(size(z))); surf(zi,'CData',ones(size(z)));
%figure; hold on; surf(z,'CData',Da); colormap jet(3); shading flat; colorbar;

return;
