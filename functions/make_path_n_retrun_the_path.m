% 폴더를 만들고, 그 폴드의 path를 받는 함수

function path = make_path_n_retrun_the_path (parentFolder,folderName)
    mkdir(parentFolder,folderName);
    path = fullfile(parentFolder,folderName);
end