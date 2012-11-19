import os
import shutil

def remove_folder(folder_name):
    """ removes all the files from the given folder;
    folder_name: name of the folder from which to remove files (it will not remove subfolders)
    remove_from_subs: if True all the subfolders will be cleared as well"""
    if os.path.exists(folder_name):
        shutil.rmtree(folder_name)
    

def remove_content(folder_name):
    """ removes all the files from the given folder;
    folder_name: name of the folder from which to remove files (it will not remove subfolders)"""    

    
    for filename in os.listdir(folder_name):
        filepath = os.path.join(folder_name, filename)
        try:
            shutil.rmtree(filepath)
        except OSError:
            os.remove(filepath)
            
            
def create_folder(folder_name):
    """ checks if given folder exists, if not it creates it"""
    
    # checks if the given folder/file exists and if not calculates the data
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
        

def file_exists(folder_name, filename):  
    """ checks if given file exists in given folder, and returns True or False"""

    filepath = os.path.join(folder_name, filename)
    exists = os.path.isfile(filepath)

    return exists      


            
            
