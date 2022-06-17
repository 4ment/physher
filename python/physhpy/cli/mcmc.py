from __future__ import annotations

import logging
import math

from .evolution import create_evolution_joint, create_evolution_parser, create_taxa

logger = logging.getLogger(__name__)


def create_mcmc_parser(subprasers):
    parser = subprasers.add_parser("mcmc", help="build a JSON file for MCMC inference")
    create_evolution_parser(parser)

    parser.add_argument(
        "--iter",
        type=int,
        default=100000,
        help="""maximum number of iterations""",
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=0,
        help="""Number of steps in stepping stone""",
    )
    parser.add_argument(
        "--stem",
        required=False,
        help="""stem for output file""",
    )
    parser.set_defaults(func=build_mcmc)
    return parser


def create_operators(json_object: dict) -> tuple[list[str], list[str]]:
    distributions = []
    parameters = []
    if isinstance(json_object, list):
        for element in json_object:
            distrs, params = create_operators(element)
            distributions.extend(distrs)
            parameters.extend(params)
    elif isinstance(json_object, dict):
        if (
            "type" in json_object
            and json_object["type"] == "parameter"
            and "!PHYSHER_FIXED" not in json_object
        ):
            lower = json_object.get("lower", -math.inf)
            upper = json_object.get("upper", math.inf)
            ref = (
                "&"
                if "dimension" not in json_object and "values" not in json_object
                else "%"
            )
            parameters.append(f"{ref}{json_object['id']}")

            if lower == 0 and upper == 1:
                distributions.append(
                    {
                        "id": "scaler",
                        "type": "operator",
                        "algorithm": "beta",
                        "x": f"{ref}{json_object['id']}",
                        "weight": 1,
                    }
                )
            else:
                distributions.append(
                    {
                        "id": "scaler",
                        "type": "operator",
                        "algorithm": "scaler",
                        "x": f"{ref}{json_object['id']}",
                        "weight": 1,
                    }
                )
                distributions.append(
                    {
                        "id": "slider",
                        "type": "operator",
                        "algorithm": "slider",
                        "x": f"{ref}{json_object['id']}",
                        "weight": 1,
                    }
                )
        elif (
            "type" in json_object
            and json_object["type"] == "simplex"
            and "!PHYSHER_FIXED" not in json_object
        ):
            distributions.append(
                {
                    "id": "dirichlet",
                    "type": "operator",
                    "algorithm": "dirichlet",
                    "x": f"${json_object['id']}",
                    "weight": 1,
                }
            )
            parameters.append(f"${json_object['id']}")
        else:
            for value in json_object.values():
                distrs, params = create_operators(value)
                distributions.extend(distrs)
                parameters.extend(params)
    return distributions, parameters


def create_multiple_mcmc(id_, mcmc, arg):
    mmcmc_dic = {
        "id": f"{id_}",
        "type": "mmcmc",
        "steps": arg.steps,
        "distribution": "beta",
        "mcmc": mcmc,
    }
    return mmcmc_dic


def create_mcmc(id_, parameters, arg):
    mcmc_dic = {
        "id": f"{id_}",
        "type": "mcmc",
        "model": "&joint",
        "length": arg.iter,
        "log": [
            {
                "id": "logger",
                "type": "logger",
                "file": f"{arg.stem}.log",
                "every": 1000,
                "models": ["&joint", "&treelikelihood"],
                "x": parameters,
            },
            {
                "id": "logger2",
                "type": "logger",
                "every": 10000,
                "models": ["&joint", "&treelikelihood"],
            },
        ],
        "operators": [],
    }
    return mcmc_dic


def build_mcmc(arg):
    if arg.clock is not None:
        arg.include_jacobian = True
    json_list = []
    taxa = create_taxa("taxa", arg)

    joint_dic = create_evolution_joint(taxa, arg)
    json_list.append(joint_dic)

    json_dic = {dic["id"]: dic for dic in json_list}

    operators, parameters = create_operators(json_list)

    mcmc_dic = create_mcmc("mcmc", parameters, arg)
    mcmc_dic["operators"] = operators

    json_dic["physher"] = []
    if arg.steps > 0:
        mmcmc_dic = create_multiple_mcmc("mmcmc", mcmc_dic, arg)
        json_dic["physher"].append(mmcmc_dic)
    else:
        json_dic["physher"].append(mcmc_dic)

    if arg.steps > 0:
        json_dic["physher"].append(
            {
                "id": "marginal",
                "type": "marginallikelihood",
                "file": f"{arg.stem}.log",
                "burnin": 100,
                "distribution": "beta",
                "steps": arg.steps,
                "algorithm": "ss",
            }
        )

    return json_dic
