def find_folders(cell_no = 1, part = 1, oscil = 0):
    """ assigns appropriate folders and filenames for different data number given
    possible combinations:
    1. (1, 1, 0), (1, 1, 1), (1, 2, 1)
    """
    
    #save_path = '/save_folder/induced/'
    save_data = '/home/maja/phdProject/analysis/swp/'
    #save_data = '/folder/saved_data/'
    save_path = '/home/maja/phdProject/data/'
    
    if cell_no == 1:
        f_dir = 'Cell 1/'
        if oscil == 0:
            save_folder = 'cell1/gap_free/'
            f_name = '25102011_0019_gap free.abf' # cell1
            intra = 1
        else:
            if part == 1:
                save_folder = 'cell1/oscil_1/'
                f_name = '25102011_0020_stim oscillo_same slice but other cell not inducing SPWs.abf'
                intra = 1
            elif part == 2:
                save_folder = 'cell1/oscil_2/'
                f_name = '25102011_0023_stim oscillo_ inducing cell.abf'
                intra = 1
        
    elif cell_no == 3:
        f_dir = 'Cell 3/'
        if oscil == 0:
            save_folder = 'cell3/gap_free/'
            f_name = '08112011_0000_gap free.abf' # cell3 (no intra)
            intra = 0
        else:
            f_name = '08112011_0017_stim oscillo.abf'
            save_folder = 'cell3/oscil/'
            intra = 1
    elif cell_no == 4:
        f_dir = 'Cell 4/'
        f_name = '08112011_0020_stim oscillo_same slice as cell 2.abf'
        save_folder = 'cell4/'
        intra = 1
    elif cell_no == 5:
        f_dir = 'Cell 5/'
        if oscil == 0:
            intra = 0
            if part == 1:
                f_name = '17112011_0000_gap free_p1.abf' # cell5 p1
                save_folder = 'cell5/part1/'
            elif part == 2:
                f_name = '17112011_0000_gap free_p2.abf'
                save_folder = 'cell5/part2/'
            elif part == 3:
                f_name = '17112011_0000_gap free_p3.abf'
                save_folder = 'cell5/part3/'
            elif part == 4:
                f_name = '17112011_0000_gap free_p4.abf'
                save_folder = 'cell5/part4/'
            elif part == 5:
                f_name = '17112011_0000_gap free_p5.abf'
                save_folder = 'cell5/part5/'
        else:
            f_name = '17112011_0002_stim oscillo.abf'
            save_folder = 'cell5/oscil/'
            intra = 1
    elif cell_no == 6:
        f_dir = 'Cell 6/'
        if oscil == 0:
            f_name = '02122011_0000_gap free stim.abf'# cell6
            save_folder = 'cell6/gap_free/'
            intra = 1
        else:
            f_name = '02122011_0001_oscillo stim.abf'
            save_folder = 'cell6/oscil/'
            intra = 1
    elif cell_no == 7:
        f_dir = 'Cell 7/'
        f_name = '02122011_0006_gap free stim.abf'
        save_folder = 'cell7/'
        intra = 1
    elif cell_no == 8:
        f_dir = 'Cell 8/'
        f_name = '08122011_0012_stim oscillo_as control take gap free recording from cell 5.abf' # cell8
        save_folder = 'cell8/'
        intra = 1
    elif cell_no == 9:
        f_dir = 'Cell 9/'        
        f_name = '08122011_0013_stim gap free_same slice as cell 4.abf'
        save_folder = 'cell9/'
        intra = 1
    elif cell_no == 10:
        f_dir = 'Cell 10/'
        f_name = '09122011_0003.abf'
        save_folder = 'cell10/'
        intra = 1
    elif cell_no == 11:
        f_dir = 'Cell 11/'
        f_name = '09122011_0004.abf'
        save_folder = 'cell11/'
        intra = 1
    read_name = save_path + f_dir + f_name
    save_folder = save_data + save_folder
    return read_name, save_folder, intra  