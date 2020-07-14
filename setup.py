from setuptools import setup

setup(
    name='sn_plotters',
    version='v1.0.0',
    description='Set of tools used to display SN pipeline results',
    url='http://github.com/lsstdesc/sn_plotters',
    author='Philippe Gris',
    author_email='philippe.gris@clermont.in2p3.fr',
    license='BSD',
    packages=['sn_plotter_simu', 'sn_plotter_fitlc','sn_plotter_metrics'],
    python_requires='>=3.5',
    zip_safe=False
)
