from distutils.core import setup

setup(
      name='meercat',
      version='0.0.1',
      packages=['meercat'],

      python_requires='>=3.7, <4',

      install_requires=['numpy','ninatool','matplotlib','scipy','cupy','cuda-cqed','open-fermion','qutip']
)