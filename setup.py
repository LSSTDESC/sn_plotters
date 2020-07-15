from setuptools import setup
# get the version here
pkg_vars  = {}

with open("version.py") as fp:
    exec(fp.read(), pkg_vars)

setup(
    name='sn_plotters',
    version= pkg_vars['__version__'],
    description='Set of tools used to display SN pipeline results',
    url='http://github.com/lsstdesc/sn_plotters',
    author='Philippe Gris',
    author_email='philippe.gris@clermont.in2p3.fr',
    license='BSD',
    packages=['sn_plotter_simu', 'sn_plotter_fitlc','sn_plotter_metrics'],
    python_requires='>=3.5',
    zip_safe=False
)
