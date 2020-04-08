# proto-covid
Python (>= 3.7) implementation of [arashbm/covid][covid] to help people
prototype models.

[covid]: https://github.com/arashbm/covid

## Getting started

Clone and install dependencies:

```bash
$ git clone https://github.com/arashbm/proto-covid.git
$ python -m venv ./venv
$ source venv/bin/activate
$ pip install --upgrade pip
$ pip install -r requirements.txt
```

Next time you need to use the repository, remember to activate the virtual
environment:

```bash
$ source venv/bin/activate
```

And when you're done working with this repository deactivate the virtual
environment:


```bash
$ deactivate
```

For an example use case scenario, check out the `example-regions.py` script.

## Model

The model is a metapopulation compartment model based on [Arenas et al.][arenas]
and revised by Eva Kisdi and others.

[arenas]: https://www.medrxiv.org/content/10.1101/2020.03.21.20040022v1
