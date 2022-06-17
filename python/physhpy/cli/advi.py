from __future__ import annotations

import logging

from .evolution import create_evolution_joint, create_evolution_parser, create_taxa

logger = logging.getLogger(__name__)


def create_variational_parser(subprasers):
    parser = subprasers.add_parser(
        "advi", help="build a JSON file for variational inference"
    )
    create_evolution_parser(parser)

    parser.add_argument(
        "--iter",
        type=int,
        default=100000,
        help="""maximum number of iterations""",
    )
    parser.add_argument(
        "-q",
        "--variational",
        nargs="*",
        default="meanfield",
        help="""variational distribution family""",
    )
    parser.add_argument(
        "--lr",
        default=0.1,
        type=float,
        help="""learning rate (default: 0.1)""",
    )
    parser.add_argument(
        "--elbo_samples",
        type=int,
        default=100,
        help="""number of samples for Monte Carlo estimate of ELBO""",
    )
    parser.add_argument(
        "--grad_samples",
        type=int,
        default=1,
        help="""number of samples for Monte Carlo estimate of gradients""",
    )
    parser.add_argument(
        "--K_grad_samples",
        type=int,
        default=1,
        help="number of samples for Monte Carlo estimate of gradients"
        " using multisample objective",
    )
    parser.add_argument(
        "--samples",
        type=int,
        default=1000,
        help="""number of samples to be drawn from the variational distribution""",
    )
    parser.add_argument(
        "--tol_rel_obj",
        type=float,
        default=0.01,
        help="""convergence tolerance on the relative norm of the objective
         (defaults: 0.001)""",
    )
    parser.add_argument(
        "--entropy",
        required=False,
        action="store_true",
        help="""use entropy instead of log Q in ELBO""",
    )
    parser.add_argument(
        "--stem",
        required=False,
        help="""stem for output file""",
    )
    parser.add_argument(
        "--convergence_every",
        type=int,
        default=100,
        help="""convergence every N iterations""",
    )
    parser.set_defaults(func=build_advi)
    return parser


def create_meanfield(var_id: str, json_object: dict) -> tuple[list[str], list[str]]:
    distributions = []
    var_parameters = []
    if isinstance(json_object, list):
        for element in json_object:
            distrs, params = create_meanfield(var_id, element)
            distributions.extend(distrs)
            var_parameters.extend(params)
    elif isinstance(json_object, dict):
        if (
            "type" in json_object
            and json_object["type"] == "parameter"
            and "!PHYSHER_FIXED" not in json_object
        ):
            if "dimension" not in json_object and "values" not in json_object:
                distr = {
                    "id": f"{var_id}.{json_object['id']}",
                    "type": "distribution",
                    "distribution": "normal",
                    "x": f"&{json_object['id']}",
                    "initialize": True,
                    "parameters": {
                        "mu": {
                            "id": f"{json_object['id']}.mu",
                            "type": "parameter",
                            "value": json_object["value"],
                        },
                        "sigma": {
                            "id": f"{json_object['id']}.sigma",
                            "type": "parameter",
                            "value": 0.1,
                            "lower": 0,
                        },
                    },
                }
                distributions.append(distr)
                var_parameters.extend(
                    (f"&{json_object['id']}.mu", f"&{json_object['id']}.sigma")
                )
            else:
                # values = np.asarray(json_object['values'])
                # if 'lower' in json_object:
                #     if json_object['lower'] != 0:
                #         values -=  json_object['lower']
                #     loc = np.log(values / np.sqrt(1 + 0.001 / values**2)).tolist()
                #     scale = np.sqrt(np.log(1 + 0.001 / values**2)).tolist()
                distr = {
                    "id": f"{var_id}.{json_object['id']}",
                    "type": "distribution",
                    "distribution": "normal",
                    "x": f"%{json_object['id']}",
                    "initialize": True,
                    "parameters": {
                        "mu": {
                            "id": f"{json_object['id']}.mu",
                            "type": "parameter",
                            "values": json_object["values"],
                            "dimension": json_object["dimension"],
                        },
                        "sigma": {
                            "id": f"{json_object['id']}.sigma",
                            "type": "parameter",
                            "values": [0.1],
                            "dimension": json_object["dimension"],
                            "lower": 0,
                        },
                    },
                }
                distributions.append(distr)
                var_parameters.extend(
                    (f"%{json_object['id']}.mu", f"%{json_object['id']}.sigma")
                )
        elif (
            "type" in json_object
            and json_object["type"] == "simplex"
            and "!PHYSHER_FIXED" not in json_object
        ):
            distr = {
                "id": f"{var_id}.{json_object['id']}",
                "type": "distribution",
                "distribution": "normal",
                "simplices": [f"%{json_object['id']}"],
                # "initialize": true,
                "parameters": {
                    "mu": {
                        "id": f"{json_object['id']}.mu",
                        "type": "parameter",
                        "values": [0],
                        "dimension": len(json_object["values"]) - 1,
                    },
                    "sigma": {
                        "id": f"{json_object['id']}.sigma",
                        "type": "parameter",
                        "values": [0.1],
                        "dimension": len(json_object["values"]) - 1,
                        "lower": 0,
                    },
                },
            }
            distributions.append(distr)
            var_parameters.extend(
                (f"%{json_object['id']}.mu", f"%{json_object['id']}.sigma")
            )
        else:
            for value in json_object.values():
                distrs, params = create_meanfield(var_id, value)
                distributions.extend(distrs)
                var_parameters.extend(params)
    return distributions, var_parameters


def create_variational_model(id_, joint, arg) -> tuple[dict, list[str]]:
    variational = {
        "id": id_,
        "type": "variational",
        "posterior": "&joint",
        "elbosamples": arg.elbo_samples,
        "distributions": [],
    }
    if arg.K_grad_samples > 1:
        variational["gradsamples"] = [arg.grad_samples, arg.K_grad_samples]
    else:
        variational["gradsamples"] = arg.grad_samples
    distributions, parameters = create_meanfield(id_, joint)
    variational["distributions"] = distributions
    return variational, parameters


def create_advi(variational, parameters, arg):
    if arg.stem:
        checkpoint = arg.stem + "-checkpoint.csv"
    else:
        checkpoint = "checkpoint.csv"
    advi_dic = {
        "id": "sg",
        "type": "optimizer",
        "algorithm": "sg",
        "update": "adam",
        "tol": arg.tol_rel_obj,
        "eta": arg.lr,
        "model": f"&{variational}",
        "checkpoint": checkpoint,
        "parameters": parameters,
        "frequency_check": arg.convergence_every,
        "max": arg.iter,
    }
    return advi_dic


def create_sampler(id_, var_id, parameters, arg):
    if arg.stem:
        file_name = arg.stem + "-samples.csv"
        tree_file_name = arg.stem + "-samples.trees"
    else:
        file_name = "samples.csv"
        tree_file_name = "samples.trees"

    parameters2 = list(filter(lambda x: "tree.ratios" != x, parameters))

    return {
        "id": id_,
        "type": "Sampler",
        "model": var_id,
        "samples": arg.samples,
        "loggers": [
            {
                "id": "logger",
                "type": "Logger",
                "file_name": file_name,
                "parameters": ["joint", "like", var_id] + parameters2,
                "delimiter": "\t",
            },
            {
                "id": "tree.logger",
                "type": "TreeLogger",
                "file_name": tree_file_name,
                "tree_model": "tree",
            },
        ],
    }


def build_advi(arg):
    if arg.clock is not None:
        arg.include_jacobian = True
    json_list = []
    taxa = create_taxa("taxa", arg)

    joint_dic = create_evolution_joint(taxa, arg)
    json_list.append(joint_dic)

    var_dic, var_parameters = create_variational_model("variational", json_list, arg)

    json_list.append(var_dic)
    json_dic = {dic["id"]: dic for dic in json_list}

    advi_dic = create_advi("variational", var_parameters, arg)
    json_list.append(advi_dic)

    json_dic["physher"] = []
    json_dic["physher"].append(advi_dic)

    parameters = []
    if arg.clock is not None:
        branch_model_id = "branchmodel"
        parameters.extend(["tree.ratios", "tree.root_height"])

        if arg.clock == "ucln":
            parameters.extend(
                (
                    f"{branch_model_id}.rates.prior.mean",
                    f"{branch_model_id}.rates.prior.scale",
                )
            )
        else:
            parameters.append(f"{branch_model_id}.rate")

        if arg.clock == "ucln":
            parameters.append(f"{branch_model_id}.rates")
    else:
        parameters = ["tree.blens"]

    if arg.coalescent is not None:
        parameters.append("coalescent.theta")
        if arg.coalescent in ("skygrid", "skyride"):
            parameters.append("gmrf.precision")
        elif arg.coalescent == "exponential":
            parameters.append("coalescent.growth")

    if arg.model == "GTR":
        parameters.extend(["substmodel.rates", "substmodel.frequencies"])
    elif arg.model == "HKY":
        parameters.extend(["substmodel.kappa", "substmodel.frequencies"])

    if arg.categories > 1:
        parameters.append("sitemodel.shape")

    # if arg.samples > 0:
    #     json_list.append(create_sampler('sampler', 'variational', parameters, arg))
    return json_dic
