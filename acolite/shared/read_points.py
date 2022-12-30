## read_points
## function to write points.txt file for plotting
## example
#label=Dumbarton Bridge
#lat=37.50472550256198
#lon=-122.11908541195334
#color=Red
#sym=o
#label_side=right
## written by Quinten Vanhellemont
## 2018-03-29
## modifications:


def read_points(file):
    points={}
    with open(file, 'r') as f:
        for il, line in enumerate(f.readlines()):
            line=line.strip()
            if len(line) == 0: continue
            if line[0] in ['#', '%', ';']: continue

            sp = line.split('=')
            if sp[0]=='label':
                cur_label = sp[1]
                if cur_label not in points:
                    points[cur_label] = {}

            if sp[1] == 'True': sp[1] = True
            if sp[1] == 'False': sp[1] = False
            if sp[1] == 'None': sp[1] = None

            points[cur_label][sp[0]]=sp[1]
    return(points)
