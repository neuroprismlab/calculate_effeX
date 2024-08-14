function quick_vis(img)

mid_z=round(size(img,3)/2); % TODO: catch error when ndims < 3
if ndims(img)==3 % 3D
    imagesc(img(:,:,mid_z));
elseif ndims(img)==4 % 4D
    imagesc(img(:,:,mid_z,1));
else
    error('Will only visualize 4D or smaller.')
end