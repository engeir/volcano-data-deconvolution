"""Main script to print info about the package."""

import vdd


def main() -> None:
    """Run the main function for the vdd package."""
    print(f"Hello, this is {__package__} at version v{vdd.__version__}!")
