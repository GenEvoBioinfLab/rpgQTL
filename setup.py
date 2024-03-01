from setuptools import setup

setup(
      name='rpgQTL',
      version='1.0.0',
      description='Regions per gene (rpg) QTL',
      url='https://github.com/gaotc200/rpgQTL',
      author='Jiahao Gao',
      author_email='jiahao.gao@yale.edu',
      license='BSD 3',
      packages=['rpgQTL'],
      zip_safe=False,
      install_requires = [
      'numpy',
      'pandas',
      'scipy',
      'torch',
      'tensorqtl']
)
