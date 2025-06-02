function [T_ap, T_pa] = pipe_distortion_correction(data, ind_ap, ind_pa)
% function [T_ap, T_pa] = pipe_distortion_correction(data)

% assume ap in first dim, pa in second volume


% 1. Logistics
op = msf_tmp_path(1);

ap_fn = fullfile(op, 'tmp_ap.nii.gz');
pa_fn = fullfile(op, 'tmp_pa.nii.gz');

I = data.imreshape();

mdm_nii_write(I(:,:,:,ind_ap), ap_fn, data.h);
mdm_nii_write(I(:,:,:,ind_pa), pa_fn, data.h);


% 2. Registration
p = elastix_p_affine(500);
p.Transform = 'BSplineTransform';
p.FinalGridSpacingInVoxels  = [ 1 1 1 ] * 20;
p_fn = elastix_p_write(p, 'p.txt');

opt.mio.coreg.clear_header = 0;
opt.do_overwrite = 1;

[~,~,T_ap] = mdm_coreg(ap_fn, pa_fn, p_fn, op, opt);
[~,~,T_pa] = mdm_coreg(pa_fn, ap_fn, p_fn, op, opt);


% 3. Logistics -> Save displacemnt fields
t_ap_fn = fullfile(op, 't_ap.txt');
elastix_p_write(T_ap.t, t_ap_fn);
msf_delete(fullfile(op, 'deformationField.nii'));
msf_system(sprintf('transformix -def all -out %s -tp %s', op, t_ap_fn));
T_ap = mdm_nii_read(fullfile(op, 'deformationField.nii'));


t_pa_fn = fullfile(op, 't_pa.txt');
elastix_p_write(T_pa.t, t_pa_fn);
msf_delete(fullfile(op, 'deformationField.nii'));
msf_system(sprintf('transformix -def all -out %s -tp %s', op, t_pa_fn));
T_pa = mdm_nii_read(fullfile(op, 'deformationField.nii'));

if (size(T_ap, 5) ~= 3), error('using an old framework? bad nii_read'); end

% 4. Edit displacement fields
T_ap = permute(T_ap, [1 2 3 5 4]) / 2;
T_pa = permute(T_pa, [1 2 3 5 4]) / 2;


