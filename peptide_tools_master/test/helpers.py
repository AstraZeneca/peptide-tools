import os


def stringify_list(lst):
    return [str(e) for e in lst]


def remove_file_list(file_list):
    for file in file_list:
        os.remove(file)


def raise_if_file_exists_list(file_list):
    for file in file_list:
        if os.path.exists(file):
            raise FileExistsError(f"File {file} already exists.")


def raise_if_file_not_exists_list(file_list):
    for file in file_list:
        if not os.path.exists(file):
            raise FileNotFoundError(f"File {file} does not exist.")
