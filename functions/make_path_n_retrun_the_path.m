% ������ �����, �� ������ path�� �޴� �Լ�

function path = make_path_n_retrun_the_path (parentFolder,folderName)
    mkdir(parentFolder,folderName);
    path = fullfile(parentFolder,folderName);
end