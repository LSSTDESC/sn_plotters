from setuptools import setup

setup(
    name='sn_plotters',
    version='0.1',
    description='Set of tools used to display SN pipeline results',
    url='http://github.com/lsstdesc/sn_plotters',
    author='Philippe Gris',
    author_email='philippe.gris@clermont.in2p3.fr',
    license='BSD',
    packages=['sn_plotter_simu', 'sn_plotter_fitlc'],
    python_requires='>=3.5',
    zip_safe=False
)
