# spooky_gui.spec
# -*- mode: python ; coding: utf-8 -*-

block_cipher = None

a = Analysis(
    ['gui/spooky_gui.py', 'core/spooky_processor.py'],
    pathex=['C:\\Users\\Bruno_Win\\Downloads\\spooky_project_20250625'],
    binaries=[],
    datas=[('README.txt', '.')],
    hiddenimports=[
        'numpy',
        'Bio._py3k',
        'Bio.SeqUtils',
        'Bio.Seq',
        'Bio.SeqRecord',
        'Bio.AlignIO',
        'Bio.Nexus',
        'Bio.Alphabet'
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        'pandas',
        'matplotlib',
        'scipy',
        'sqlalchemy'
    ],
    noarchive=False,
    cipher=block_cipher,
    #target_arch=None,
    #upx=True,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='spooky_gui',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    target_arch=None,
)