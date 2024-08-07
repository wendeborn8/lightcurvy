from setuptools import setup

with open('requirements.txt') as f:
    install_reqs = f.read().splitlines()

setup(
      name = 'Lightcurvy',
      version = '0.2',
      description = 'Scripts to fetch astronomical light curves from a variety of sources.',
      author = 'John Carlos Wendeborn',
      author_email = 'jwendeborn@gmail.com',
      install_requires = install_reqs,
      license = 'MIT',
      url = 'https://github.com/wendeborn8/lightcurvy',
      )
