function calc_fc_hcpU100(type)

    addpath(genpath('/n02dat01/users/yfwang/software/MatlabToolbox/NIFTI'));
    addpath(genpath('/n02dat01/users/yfwang/software/MatlabToolbox/spm12/toolbox/brant/brant_postprocess/brant_STAT'));
    addpath(genpath('/share/soft/brant'));
    datapath='/n04dat01/atlas_group/yfwang/HCP_U100/';
    resultpath='/n02dat01/users/ypwang/Gradient/GeneticGradient/Data/Sigsample2CereFC_lesionnetwork/Data_HCP_U100/';
    
    ref = load_nii(strcat(resultpath, 'ref_float_2mm.nii.gz'));
    ref_nii = ref.img;

    % read subject list
    subs = textread(strcat(datapath, 'HCP_U100_list.txt'), '%s');

    % generate GM coord file
    targetcoordfile = strcat(resultpath, 'cerebellum_2mm_mask.txt');
    if ~exist(targetcoordfile) 
        target = load_nii(strcat(resultpath, 'Cerebellum_2mm_mask_17853.nii'));
        target.img(isnan(target.img))=0;
        [nx,ny,nz] = size(target.img);
        fid = fopen(targetcoordfile, 'w');
        for z = 1:nz
            [x y] = find(target.img(:,:,z) == 1);
            for j = 1:numel(x)
                fprintf(fid,'%d %d %d\r\n',x(j),y(j),z);
            end
        end
        fclose(fid);
    end
    targetcoordinates = int16(importdata(targetcoordfile));
    
    % calc FC
    if ~exist(strcat(resultpath, type, '_z.mat'))
        coordinates = load(strcat(resultpath, type, '.txt'));
        
        corr_r = zeros([length(subs), length(targetcoordinates)]); % 100*17853
        corr_z = zeros([length(subs), length(targetcoordinates)]); % 100*17853

        for j=1:length(subs)
            disp(j)
            namelist = ['s6fdGSRrfMRI_REST1_LR_hp2000_clean'; 's6fdGSRrfMRI_REST1_RL_hp2000_clean'; 's6fdGSRrfMRI_REST2_LR_hp2000_clean'; 's6fdGSRrfMRI_REST2_RL_hp2000_clean'];
            len = 4;
            for i = 1:len
                fmri = load_nii(strcat(datapath, subs{j}, '/', namelist(i,:), '.nii.gz'));

                gm = zeros([length(targetcoordinates), size(fmri.img,4)]); % 17853*1200
                for k=1:length(targetcoordinates)
                    x=targetcoordinates(k,1);
                    y=targetcoordinates(k,2);
                    z=targetcoordinates(k,3);
                    gm(k,:) = transpose(squeeze(fmri.img(x,y,z,:)));
                end

                ts = zeros([length(coordinates), size(fmri.img,4)]); % ~150(vertices)*1200
                for k=1:length(coordinates)
                    x=coordinates(k,1);
                    y=coordinates(k,2);
                    z=coordinates(k,3);
                    ts(k,:)=transpose(squeeze(fmri.img(x,y,z,:)));
                end

                meants=mean(ts); % 1*1200
                tmp_r = corr(meants', gm');
                corr_r(j,:) = corr_r(j,:) + tmp_r;

                clear fmri;
                clear gm;
                clear tmp_r;
            end
            corr_r(j,:) = corr_r(j,:)./len;
        end
        corr_z = (log((1+corr_r) ./ (1-corr_r))) / 2;
        
        save(strcat(resultpath, type, '_r.mat'), 'corr_r');
        save(strcat(resultpath, type, '_z.mat'), 'corr_z');
    end

    % ttest
    cere_coord = targetcoordinates;
    cere_num = size(cere_coord,1); % 17853

    corr_z = load(strcat(resultpath, type, '_z.mat'));
    corr_z = corr_z.corr_z;
    p = zeros([1,cere_num]);
    t = zeros([1,cere_num]);
    for j=1:cere_num
        [h,p(j),ci,st]=ttest(corr_z(:,j));
        t(j) = st.tstat;
    end
    save(strcat(resultpath, type, '_ttest.mat'), 'h','p','ci','t');
    
    timg = zeros(size(ref_nii));
    for j=1:cere_num
        x = cere_coord(j,1);
        y = cere_coord(j,2);
        z = cere_coord(j,3);
        timg(x,y,z) = p(j);
    end
    ref.img = timg;
    save_nii(ref, strcat(resultpath, type, '_p.nii.gz'));

    timg = zeros(size(ref_nii));
    for j=1:cere_num
        x = cere_coord(j,1);
        y = cere_coord(j,2);
        z = cere_coord(j,3);
        timg(x,y,z) = t(j);
    end
    ref.img = timg;
    save_nii(ref, strcat(resultpath, type, '_t.nii.gz'));

    % fdr/fwe
    thres = [0.001, 0.01, 0.05];

    data = load(strcat(resultpath, type, '_ttest.mat'));
    p = data.p;
    t = data.t;
    
    % fdr
    for thr=thres
        timg = zeros(size(ref_nii));
        [pfdr, b] = brant_MulCC(p',thr,'fdrID');
        for j=1:cere_num
            x = cere_coord(j,1);
            y = cere_coord(j,2);
            z = cere_coord(j,3);
            if p(j) < pfdr
                timg(x,y,z) = t(j);
            end
        end
        ref.img = timg;
        save_nii(ref, strcat(resultpath, type, '_t_fdr_', num2str(thr), '.nii.gz'));
    end

    % fwe
    for thr=thres
        timg = zeros(size(ref_nii));
        [pfwe, b] = brant_MulCC(p',thr,'bonf');
        for j=1:cere_num
            x = cere_coord(j,1);
            y = cere_coord(j,2);
            z = cere_coord(j,3);
            if p(j) < pfwe
                timg(x,y,z) = t(j);
            end
        end
        ref.img = timg;
        save_nii(ref, strcat(resultpath, type, '_t_fwe_', num2str(thr), '.nii.gz'));
    end
    disp('done')
