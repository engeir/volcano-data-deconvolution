# Accompanying code to 'Non-parametric estimation of temperature response to volcanic forcing'

[![DOI](https://zenodo.org/badge/742292778.svg)](https://zenodo.org/badge/latestdoi/742292778)

<sup>Latest version: v2.0.1</sup> <!-- x-release-please-version -->

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
environment into `./.venv/`. If you do not wish to use [mise], [uv] will by default also
create a virtual environment with an available python runtime.

```bash
git clone git@github.com:engeir/volcano-data-deconvolution.git
cd volcano-data-deconvolution || exit
uv sync
```

> [!NOTE]
>
> The repository and project is named volcano-data-deconvolution, but the package name
> is `vdd`, just to make it shorter when importing.

After the `uv sync` command succeeds, you can check that everything is working by
running

<!-- x-release-please-start-version -->

```console
$ uv run vdd
Hello, this is vdd at version v2.0.1!
```

<!-- x-release-please-end -->

## Usage

The main programs, i.e. those that create plots used in the accompanying paper, are run
as entry point scripts of the package. A list of all entry points can be found by
running the command `uv run`, or for a curated list of this project's entry points, run

```bash
uv run | grep vdd
```

All other useful scripts are located in the [`src/vdd/plotting/`](./src/vdd/plotting/)
directory. Any of those scripts can be run as

```bash
uv run python ./src/vdd/plotting/some_plotting_script.py
```

> [!NOTE]
>
> And again, if you use mise, you can do
>
> ```bash
> mise run plot
> ```
>
> and it will magically just work.

[mise]: https://mise.jdx.dev/
[uv]: https://docs.astral.sh/uv/
