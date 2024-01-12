# Deconvolution of volcanic eruption time series

<sup>Latest version: v0.1.2</sup> <!-- x-release-please-version -->

> [!WARNING]
>
> This README reflects the changes made to the main branch. For the most up to date
> documentation about the version you are using, see the README at the relevant tag.

This repository contains analysis scripts that analyse output from volcanic eruption
simulations in climate models, using a deconvolution algorithm commonly applied to time
series with strong intermittency.

## Install

To install the project you must clone the repository. If you have [mise] installed
(recommended), this will install and set up the correct python version as a virtual
environment into `./.venv/`. If you do not wish to use [mise], [poetry] will by default
also create a virtual environment with an available python runtime.

```bash
git clone git@github.com:engeir/volcano-data-deconvolution.git
cd volcano-data-deconvolution
# poetry can be installed from pipx, as `pipx install poetry`. See https://python-poetry.org/docs/#installation
poetry install
```

> [!NOTE]
>
> The repository and project is named volcano-data-deconvolution, but the package name
> is `vdd`, just to make it shorter when importing.

After the `poetry install` command succeeds, you can check that everything is working by
running

<!-- x-release-please-start-version -->

```console
$ poetry run vdd
Hello, this is vdd at version v0.1.2!
```

<!-- x-release-please-end -->

## Usage

So far, only the welcome message at `poetry run vdd` is implemented. More will come!

[poetry]: https://python-poetry.org
[mise]: https://mise.jdx.dev/
