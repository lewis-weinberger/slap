# slap

> **S**imulating **L**yman-**A**lpha emitting galaxy **P**opulations

A Python implementation of the LAE modelling described in [Weinberger et al. 2019](https://doi.org/10.1093/mnras/stz481).

## Installation

**slap** has the following dependencies:

- [Python 3](https://www.python.org/downloads/)
- [NumPy](https://www.numpy.org/)
- [SciPy](https://www.scipy.org/)
- [hmf](https://github.com/steven-murray/hmf)
- [h5py](https://www.h5py.org/) (optional)
- [matplotlib](https://matplotlib.org/) (optional)

The optional h5py and matplotlib libraries are used in the provided example, [example.py](./example.py). A guide on how to set up a [virtualenv](https://virtualenv.pypa.io/en/latest/) and install the above libraries can be found in [SETUP](./SETUP.md).

The **slap** code can be downloaded using git:

```
$ git clone https://github.com/lewis-weinberger/slap
```

## Usage

**slap** provides two classes: `LAEModel` (defined in [slap.py](./slap.py)) and `Catalogue` (defined in [catalogue.py](./catalogue.py)). The `LAEModel` class can be used to generate an LAE population from a dark matter halo population. The `Catalogue` class is a useful data-structure for storing the generated LAE population in a catalogue.

A full example detailing the use of **slap** can be found at [EXAMPLE](./EXAMPLE.md).

## License

[MIT License](./LICENSE)
