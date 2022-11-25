function B2_delete_scans(run_dir, filter, nr_scans)

f=spm_select('List', run_dir, filter);

file=[run_dir filesep f];
%spm_file_split(file);

N=load_untouch_nii(file);
N.img(:,:,:,1:nr_scans)=[];
N.hdr.dime.dim(5)=N.hdr.dime.dim(5)-nr_scans;
save_untouch_nii(N,[run_dir filesep 'x' num2str(nr_scans) f]);
%
% delete(file);
%
% files=spm_select('List', run_dir, filter);
% numVols = size(files,1);
% files=[repmat([run_dir filesep], numVols, 1) files];
%
% AllVolume=[];
% for i=1:nr_scans
%     delete(files(i,:))
% end
%
% files=spm_select('List', run_dir, filter);
%
% matlabbatch{1}.spm.util.cat.vols = files';
% matlabbatch{1}.spm.util.cat.name = [run_dir filesep 'a' f];
% matlabbatch{1}.spm.util.cat.dtype = 4;
%
% %change_3D_to_4D(run_dir, filter, f);
%
% delete(files)
%
% % for i=(nr_scans+1):numVols
% %     delete(files(i,:));
% % end