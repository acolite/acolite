## def read_list
## read in text file as list
## written by Quinten Vanhellemont, RBINS
## 2022-03-28
## Last modifications:

def read_list(file):
    data_list = []
    with open(file, 'r', encoding='utf-8') as f:
        for il, line in enumerate(f.readlines()):
            line = line.strip()
            if line[0] in ['#', ';']: continue
            data_list.append(line)
    return(data_list)
