function[image_cell,pixelSize,window_width]=read_dicom_series(input_folder,save_folder)
% Function reading DICOM files from an input folder of the type "I(number)"
% and returning as outputs all the series of images present in that DICOM
% set saved in another folder. In that folder, there is one folder per
% series found. It generates a cell containing a 3D stack of images per
% series generated and a matrix with all the pixel sizes of the images
% found. The window width is also extracted for obtaining the optimal
% window width for representing CT images.

addpath(input_folder);
st1=mkdir(input_folder);
addpath(save_folder);
st2=mkdir(save_folder);
if st1==0 % Looks if the input folder where the images exists or does not exist. 
    % In case it does not exist, it creates the folder and adds the path to
    % the directory.
    mkdir(input_folder);
    addpath(input_folder);
end

if st2==0 % Looks if the input folder where the images exists or does not exist. 
    % In case it does not exist, it creates the folder and adds the path to
    % the directory.
    mkdir(save_folder);
    addpath(save_folder);
end



%%%% Obtain the number of DICOM files in the input folder %%%%
% We assume that all the DICOM files start with 'I' as first character
list=dir(input_folder); % Obtain all the files features in the input folder
no_dicom=0; % Counter of not DICOM files
for i=1:length(list)
    if strcmp(list(i).name(1),'I')==0
        no_dicom=no_dicom+1;
    end
end
n=length(list)-no_dicom; % Number of DICOM files in the input folder



if n>0 
    % If any DICOM file has been found in the input folder (starting by 'I')
    
    cont=1; % It will be used as index to locate different features of the set of DICOM images. Only counts  DICOM files
    patient={}; % Cell of length n containing the patient IDs of all files
    patient_norep={}; % Cell containing the patient IDs found in all DICOM headers without repetition
    studies={}; % Cell of length n containing all the studies UIDs found in the DICOM set of images
    studies_norep={}; % Cell containg all the studies UIDs found without repetition
    series={}; % Cell of length n containg all the series UIDs found in the DICOM set of images
    series_norep={}; % Cell containg all the series UIDs found without repetition
    thick=[]; % Vector of length n containing all the slice thickness
    pixel_size=[]; % Vector containing the pixel size of the DICOM images
    modality={}; % Cell of length n containing the modalities of image in the set
    window_width=[]; % Matrix of n x 2 size containing the window widths of each image
    rescale_slope=zeros; % Vector containing the rescale slope of each CT image to obtain the intensity values in Hounsfield Units
    rescale_intercept=zeros;% Vector containing the rescale intercept of each CT image to obtain the intensity values in Hounsfield Units
    patient_position=[];% Keep the patient position to later on know the indexes of the slices in case the series is a 3D volume
    cont_full=0; % Counter for full files
    
    while cont_full<n % Go through all images in the DICOM set. Until all the full DICOM files have not been read, it will not stop
        if(exist(strcat('I',num2str(cont-1)))~=0) % Looks if an image in DICOM format in the input folder exists or not
            % As the name of the DICOM files goes from I0 and on, we can look for
            % those files by using the index stated before, converting it to a
            % string and concatenating it to I
            info=dicominfo(strcat('I',num2str(cont-1))); % Obtain the header of each DICOM file
            info.cont=cont; % New field in the DICOM header of each image in order to identify it with a counter
            modality{cont}=info.Modality;
            if strcmp(modality{cont},'CT')==1 % If the modality is CT, express the intensity values in Hounsfield Units
                window_width(:,cont)=info.WindowWidth*info.RescaleSlope+info.RescaleIntercept;
            else
                window_width(:,cont)=info.WindowWidth;
            end
            patient{cont}=info.PatientID;
            studies{cont}=info.StudyInstanceUID;
            series{cont}=info.SeriesInstanceUID;
            pixel_size(cont)=info.PixelSpacing(1);
            thick(cont)=info.SliceThickness;
            rescale_slope(cont)=info.RescaleSlope;
            rescale_intercept(cont)=info.RescaleIntercept;
            patient_position(:,cont)=info.ImagePositionPatient;
            cont=cont+1;
            cont_full=cont_full+1;
        else
            % In case some intermediate DICOM file has been deleted
            % but there are still DICOM files unread in the input folder, it
            % goes number by number increasing the counter until the next DICOM
            % file is found. The image is information is left as zero.
            % We have to specify a special counter for empty files: cont_empty
            modality{cont}='No';
            patient{cont}='No';
            studies{cont}='No';
            series{cont}='No';
            pixel_size(cont)=0;
            thick(cont)=-1; % Leave it as negative as there could be 2D acquisitions with no thickness
            window_width(:,cont)=zeros(1,2);
            cont=cont+1;
        end
    end
    cont=cont-1; % Compensate the extra cont that has been added after the while loop
    % Obtain the absolute value of the patient positions
    patient_position=abs(patient_position);
    % Obtain with the command "unique" the different patients, studies and
    % series obtained
    patient_norep=unique(patient);
    studies_norep=unique(studies);
    series_norep=unique(series);

    
    % Now we look for each patient in the cell without repetitions and look
    % for the indexes in the counter where the patient appears. Those
    % indexes will indicate which images belong to every patient
    
    patient_index=cell(length(patient_norep)); % Cell with the length of patient_norep containing for each patient the 
    % vectors with the indexes where each patient is located in the DICOM
    % set of images
    study_index={}; % Cell with the indexes of the DICOM set of images of the studies in studies_norep being divided in patients and inside each
    % patient, the study that has been conducted
    series_index={}; % Cell with the indexes of the DICOM set of images corresponding to the series in series_norep
    % being divided in patient, the study that has been conducted and the
    % series inside each study
    aux_patient=[]; % Auxiliary vector for working with patient cell
    aux_studies=[]; % Auxiliary vector for working with studies cell
    aux_series=[]; % Auxiliary vector for working with series
    study_index_perPatient={}; % Cell with the studies organized per patient
    series_index_perStudy={}; % Cell with the series organized per study
    series_index_perPatient={}; % Cell with the series organized per study and patient
    count_image_series=1; % Counter for the image inside each series
    image=[]; % Image read with dicomread from the input folder and assigned to a certain DICOM file in the saved folder
    matrix3D=[]; % 3D matrix with slices without ordering
    image3D=[]; % Matrix keeping the values of the read DICOM images in 3D and with those slices ordered
    image_cell={}; % Cell with image information organized in 3D matrices
    position_submatrix=zeros(3,n); % Submatrix of the already obtained patient positions to get the index of the slice in each series image
    aux_position=ones; % Vector with each row of the matrix above
    correct_index=0; % Correct index for the order of the slices in each series
    
    for i=1:length(patient_norep)
           if strcmp(patient_norep{i},'No')==0 % Do not take into account empty patients
                mkdir(strcat(save_folder,'\Patient ',num2str(i))); % Create a folder per patient inside the save folder
                for j=1:cont
                    if strcmp(patient_norep{i},patient{j})==1
                        aux_patient=[aux_patient j];
                    end
                end
                patient_index{i}=aux_patient;
                aux_patient=[];
                if length(patient_norep)<=length(studies_norep) % Look if we have less patients than studies or as patients as studies
                    for l=1:length(studies_norep)
                        if strcmp(studies_norep{l},'No')==0 % Do not take into account empty studies
                            for k=patient_index{i} % Now look for the studies of each patient
                                if strcmp(studies_norep{l},studies{k})==1
                                    aux_studies=[aux_studies k];
                                end
                            end
                            if isempty(aux_studies)==0 % Look if the studies assessed belong to the patient in the iteration. If the vector is empty, no studies belong to that patient, so no folders are created
                                mkdir(strcat(save_folder,'\Patient ',num2str(i),'\Study ',num2str(l))); % Create a folder per study inside each patient folder
                                study_index_perPatient{l}=aux_studies; % Cell keeping the indexes of each study per patient
                            end
                            aux_studies=[];
                            if length(series_norep)>=length(studies_norep)
                                for m=1:length(series_norep)
                                    if strcmp(series_norep{m},'No')==0
                                        for o=study_index_perPatient{l}
                                            if strcmp(series_norep{m},series{o})==1
                                                aux_series=[aux_series o];
                                            end
                                        end
                                        if isempty(aux_series)==0 % Look if the studies assessed belong to the patient in the iteration. If the vector is empty, no series belong to that study, so no folders are created
                                            mkdir(strcat(save_folder,'\Patient ',num2str(i),'\Study ',num2str(l),'\Series ',num2str(m))); % Create a folder per study inside each patient folder
                                            for p=aux_series % Save our images in the folder
                                                % In here, we relate each image in the
                                                % input folder to the indexes obtained
                                                % knowing that the image I(index-1)
                                                % corresponds to the index obtained in here
                                                % and saved in the vector aux_series
                                                image=dicomread(strcat('I',num2str(p-1)));
                                                if strcmp(modality{p},'CT')==1
                                                    image=double(image); % Convert the image to double in order to convert it to Hounsfield Units
                                                    image=rescale_slope(p)*image+rescale_intercept(p);
                                                end
                                                position_submatrix(:,p)=patient_position(:,p);
                                                matrix3D(1:size(image,1),1:size(image,2),count_image_series)=image;
                                                count_image_series=count_image_series+1;
                                            end
                                            for q=1:size(patient_position,1)
                                                aux_position=position_submatrix(q,aux_series)-position_submatrix(q,aux_series(1));
                                                if norm(aux_position)~=0 % We know that the Image Position Patient changes in one direction and in that direction, they change in the same order
                                                    position_submatrix(q,aux_series)=sort((position_submatrix(q,aux_series)));
                                                    for r=1:length(aux_series) % Obtain the slices in the correct order
                                                        correct_index=find(position_submatrix(q,aux_series)==patient_position(q,aux_series(r)));
                                                        image=matrix3D(:,:,r);
                                                        image3D(1:size(matrix3D,1),1:size(matrix3D,2),correct_index)=image;
                                                        dicomwrite(image,strcat('Series ',num2str(m),' Slice ',num2str(correct_index)));
                                                        movefile(strcat('Series ',num2str(m),' Slice ',num2str(correct_index)),strcat(save_folder,'\Patient ',num2str(i),'\Study ',num2str(l),'\Series ',num2str(m)));
                                                    end
                                                end
                                            end
                                            aux_position=[];
                                        end
                                        % Obtain a matrix with the pixel size for each
                                        % series and a cell with the 3D images per series
                                        pixelSize(m,1:2)=zeros(1,2)+pixel_size(aux_series(1)); % Assume the pixel size is constant in the 3D image
                                        pixelSize(m,3)=thick(aux_series(1)); % Assume the slice thickness is the third dimension of the pixel size
                                        % Assume the slice thickness is constant through
                                        % the 3D image
                                        image_cell{m}=image3D;
                                        image3D=[];
                                        count_image_series=1;
                                        series_index_perStudy{m}=aux_series;
                                        aux_series=[];
                                    end
                                end
                            else
                                disp('There are more studies than series. Some series are lacking or some series are being unread');
                            end
                            series_index_perPatient{l}=series_index_perStudy;
                        end
                    end
                    study_index{i}=study_index_perPatient;
                    series_index{i}=series_index_perPatient;
                else
                    disp('There are more patients than studies. Some patients lacks studies or some studies are being unread');
                end
           end
    end
    
    
else
    fprintf('No DICOM files found in %s',input_folder);
end


end
    