function S = csf_spatiotemp_kelly( rho, nu )
% rho - spatial frequency in cpd
% nu - temporal frequency in Hz

S = csf_spatiovel_kelly( rho, nu./rho );

end