from setuptools import setup

setup(
      name = 'Lightcurvy',
      version = '0.2',
      description = 'Scripts to fetch astronomical light curves from a variety of sources.',
      author = 'John Carlos Wendeborn',
      author_email = 'jwendeborn@gmail.com',
      install_requires = ['numpy', 'pandas', 'matplotlib', 'astroquery', 'astropy', 'tglc', 'eleanor', 'lightkurve', 'tess-point', 'scikit-learn', 'requests'],
      license = 'MIT',
      url = 'https://github.com/wendeborn8/lightcurvy',
      
      )