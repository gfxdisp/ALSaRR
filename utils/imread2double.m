function img = imread2double( file_name )

img = imread( file_name );

if( isa( img, 'uint8' ) )
    img = double(img) / 255;
elseif( isa( img, 'uint16' ) )
    img = double(img) / (2^16-1);
elseif( isa( img, 'logical' ) )
    img = double(img);    
else
    error( 'unknown data format' );
end

assert( all( img(:)<=1 ) );

end