import argparse
import json

from physhpy.cli.advi import create_variational_parser
from physhpy.cli.mcmc import create_mcmc_parser
from physhpy.cli.optimize import create_maximize_parser


def remove_unwanted(obj):
    if isinstance(obj, list):
        for element in obj:
            remove_unwanted(element)
    elif isinstance(obj, dict):
        for key in list(obj.keys()).copy():
            if not key.startswith("!PHYSHER"):
                remove_unwanted(obj[key])
            else:
                del obj[key]


def main():
    parser = argparse.ArgumentParser(
        prog="physhpy-cli",
        description="Command line interface for creating JSON file for physher",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    subprasers = parser.add_subparsers()

    create_variational_parser(subprasers)

    create_maximize_parser(subprasers)

    create_mcmc_parser(subprasers)

    arg = parser.parse_args()
    json_dic = arg.func(arg)

    remove_unwanted(json_dic)

    print(json.dumps(json_dic, indent=2))
