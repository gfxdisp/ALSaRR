function [rho_x, rho_y] = create_rho_2D( im_size, nf_rows, nf_cols )
% Create frequency ceoffs
%
% Useful for constructing Fourier-domain filters based on OTF or CSF data.
%

half_size = floor(im_size/2);
odd = mod( im_size, 2 );
freq_step = [nf_rows nf_cols]./half_size;

if( odd(2) )
    xx = [ linspace( 0, nf_cols, half_size(2)+1 ) linspace( -nf_cols, -freq_step(2), half_size(2) ) ];
else
    xx = [ linspace( 0, nf_cols-freq_step(2), half_size(2) ) linspace( -nf_cols, -freq_step(2), half_size(2) ) ];
end

if( odd(1) )
    yy = [ linspace( 0, nf_rows, half_size(1)+1 ) linspace( -nf_rows, -freq_step(1), half_size(1) ) ];
else
    yy = [ linspace( 0, nf_rows-freq_step(1), half_size(1) ) linspace( -nf_rows, -freq_step(1), half_size(1) ) ];
end

[rho_x, rho_y] = meshgrid( xx, yy );

end
