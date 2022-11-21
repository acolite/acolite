## def metadata_parse
## parse (old?) IKONOS metadata, not perfect, but good enough?
## written by Quinten Vanhellemont, RBINS
## 2022-09-15
## modifications: 2022-09-15 (QV) made all metadata tags lists, in case there are multiple components

def metadata_parse(metafile):
    keep = False
    main_tag, sec_tag = None, None
    cur_data = {}
    metadata = {}

    with open(metafile, 'r', encoding='utf-8') as f:
        for il, line in enumerate(f.readlines()):
            if line[0] == '=':
                next_section = True
            else:
                next_section = False

            if (next_section):
                if (len(cur_data) == 0): continue
                for k in cur_data:
                    metadata[k] = cur_data[k]
                cur_data = {}
                main_tag, sec_tag = None, None
                #stop
            if line[0] in ['=', '-']: continue

            white_len = len(line) - len(line.strip()) -1
            if line[0] == ' ':
                level = 1
            else:
                level = 0

            line = line.strip()
            if len(line) == 0: continue

            ## skip first part of data
            if (keep is False) & ('Product Order Metadata' in line): keep = True
            if not keep: continue

            #print(level, line)
            sp = None
            if ':' in line:
                spos = line.find(':')
                sp = [s.strip() for s in [line[0:spos], line[spos+1:]]]
                #sp = [s.strip() for s in line.split(':')]

            if sp is None:
                ## if at 0 level white space and no main tag is set, set it here
                if (level == 0) & (main_tag is None):
                    main_tag = line
                ## else set secondary tag as part of main tag
                else:
                    sec_tag = line
                #print(main_tag, sec_tag)
                continue

            ## if back at 0 whitespace remove header
            if level == 0: sec_tag = None

            if sec_tag is None:
                if main_tag not in cur_data: cur_data[main_tag] = {}
                tar = cur_data[main_tag]
            #elif (len(header) > level):
            #    if (header[level] not in cur_data[main_tag]):
            #        cur_data[main_tag][header[level]] = {}
            #    tar = cur_data[main_tag][header[level]]
            else:
                if sec_tag not in cur_data[main_tag]:
                    cur_data[main_tag][sec_tag] = {}
                tar = cur_data[main_tag][sec_tag]

            if sp[0] not in tar:
                if len(sp) == 2:
                    tar[sp[0]] = [sp[1]]
                elif len(sp) > 2:
                    tar[sp[0]] = sp[1:]
            else:
                if type(tar[sp[0]]) is not list:
                    tar[sp[0]] = [tar[sp[0]]]
                #if len(sp) == 2:
                #    cur_data[tag][sp[0]].append(sp[1])
                #elif len(sp) > 2:
                #    cur_data[tag][sp[0]] += sp[1:]
                tar[sp[0]] += sp[1:]
            #print(cur_data)
    return(metadata)
