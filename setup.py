import sys
from cx_Freeze import setup, Executable

include_files = []

build_exe_options = {'packages': [''], 'excludes': []}

base = None

# if sys.platform == 'win32':
#     base = 'Win32GUI'

exe = Executable(script = 'Cuco.py', base = base, icon = None)

setup(name = 'Cuco_executável',
      version = 'Beta_0.1',
      description = 'Cálculo de Configuração para Recarga de Reator PWR',
      options = {'build_exe': {'include_files': include_files}},
      executables = [exe])