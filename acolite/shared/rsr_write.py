## def rsr_write
## exports rsr file
## written by Dirk Aurin, NASA Goddard Space Flight Center
## 2021-07-20

def rsr_write(file, header, sensor, rsr):

    if file is not None:
        with open(file, 'w', encoding='utf-8') as f:            
            for line in header:
                f.write(line + '\n')
            f.close

        with open(file, 'a', encoding='utf-8') as f: 
            for band in rsr:
                f.write(f";; {sensor} Band: {band}" + '\n')
                for i, wave in enumerate(rsr[band]['wave']):
                    # f.write(f"{wave}\t{rsr[band]['response'][i]}\n")
                    f.write('{:.4f}\t{}\n'.format(wave,rsr[band]['response'][i]))
            f.close