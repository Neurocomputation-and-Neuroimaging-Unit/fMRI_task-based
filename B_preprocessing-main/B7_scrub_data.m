function B7_scrub_data(data_dir, filter, outliers, prefix)
%replaces outlier scans with average of neaerest 'good' neighbors

f = spm_select('List', data_dir, filter);
rest_file = fullfile(data_dir, f);

if sum(outliers)~=0
    display('load data ... ')
    [data, header]=rest_ReadNiftiImage([rest_file]);
    header.fname=[data_dir filesep prefix f];
    nearest_neighbor_pre=[];
    
    fprintf('scrubbing.')
    for i=1:length(outliers)
        fprintf('.')
        nearest_neighbor_post=[];
        if outliers(i)==0
            nearest_neighbor_pre=i;
        elseif outliers(i)==1
            %find neighbor post
            if i<=size(outliers)
                for j=(i+1):(length(outliers)-i)
                    if outliers(j)==0 && isempty(nearest_neighbor_post)
                        nearest_neighbor_post=j;
                    end
                end
            end
            if isempty(nearest_neighbor_post)
                data(:,:,:,i)=data(:,:,:,nearest_neighbor_pre);
            elseif isempty(nearest_neighbor_pre)
                data(:,:,:,i)=data(:,:,:,nearest_neighbor_post);
            else
                data(:,:,:,i)=(data(:,:,:,nearest_neighbor_pre)+data(:,:,:,nearest_neighbor_post))/2;
            end
        end
    end
    fprintf('\nsave data ... ')
    rest_Write4DNIfTI(data, header, header.fname);
else
    fprintf('no scrubbing --> copy data ... ')
    copyfile(rest_file,[data_dir filesep prefix f]);
end
display('done')

