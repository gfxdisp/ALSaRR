function S = stcsf_cam_daly(rho, nu, L_bkg,sccsf)
% rho - spatial frequency in cpd (MxN)
% nu - temporal frequency in Hz (MxN)
% L_bkg - luminance in cd/m^2 (MxN)

S_st = csf_spatiotemp_daly( rho(:), nu(:) )./(csf_spatiotemp_daly( rho(:), zeros(size(rho(:))) ) + 1e-5);
LMS_d65 = xyz2lms2006( whitepoint( 'd65' ) );
L_bkg = ones(size(rho)) .* L_bkg;
sigma = 1/5;%-1./rho(:);
A_cm = pi*(sigma).^2;
S_sp = sccsf.sensitivity_coldir( rho(:), L_bkg(:) * LMS_d65, 1, A_cm );
S = reshape( S_sp .* S_st, size(L_bkg) );
end

