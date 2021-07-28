## .spec file for generating PyInstaller binaries
##
## run with python -m PyInstaller acolite.spec
##
## run with python -m PyInstaller --noconfirm acolite.spec
## --noconfirm will overwrite existing build/dist directories

# -*- mode: python -*-
block_cipher = None

datas = []
hiddenimports = []
excludes = []

## not found
## hiddenimports+=['pandas._libs.tslibs.timedeltas']
## hiddenimports+=['pandas._libs.tslibs.np_datetime','pandas._libs.tslibs.nattype','pandas._libs.skiplist']

## hidden imports
hiddenimports+=['scipy._lib.messagestream']
hiddenimports+=['cftime']
hiddenimports+=['pywt._extensions._cwt']
hiddenimports+=['scipy.spatial.transform._rotation_groups', 'cmath']
hiddenimports+=['pyhdf.six']

a = Analysis(['launch_acolite.py'],
             binaries=[],
             datas=datas,
             hiddenimports=hiddenimports,
             hookspath=[],
             runtime_hooks=[],
             excludes=excludes,
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='acolite',
          debug=False,
          strip=False,
          upx=True,
          console=True)

coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='acolite')
