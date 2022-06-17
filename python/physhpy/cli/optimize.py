from __future__ import annotations

from .evolution import (
    create_evolution_joint,
    create_evolution_parser,
    create_taxa,
    create_tree_likelihood,
)


def create_maximize_parser(subprasers):
    parser = subprasers.add_parser(
        "maximize", help="build a JSON file for maximizing things"
    )
    create_evolution_parser(parser)

    parser.add_argument(
        "--max_iter",
        type=int,
        default=10000,
        help="""maximal number of iterations per optimization step (default: 10000)""",
    )
    parser.add_argument(
        "--stem",
        required=True,
        help="""stem for output files""",
    )
    parser.add_argument(
        "--map",
        action="store_true",
        help="""perform maximum a posteriori""",
    )
    parser.set_defaults(func=build_optimizer)
    return parser


def create_meta_optimizer(joint, arg):
    meta = {
        "id": "metaopt",
        "type": "optimizer",
        "algorithm": "meta",
        "precision": 0.001,
        "max": arg.max_iter,
        "model": f"&{joint}",
        "list": [],
    }
    list_opt = []

    if arg.clock is None:
        list_opt.append(
            {
                "id": "optbl",
                "type": "optimizer",
                "algorithm": "serial",
                "model": f"&{joint}",
                "_parameters": ["%tree.distances"],
                "treelikelihood": f"&{joint}",
                "precision": 0.00001,
            }
        )

    subst_parameters = None
    if arg.model == "GTR":
        subst_parameters = ["$frequencies", "$rates"]
    elif arg.model == "HKY":
        subst_parameters = ["$frequencies", "&kappa"]
    if subst_parameters is not None:
        list_opt.append(
            {
                "id": "opt.subst",
                "type": "optimizer",
                "algorithm": "serial",
                "model": f"&{joint}",
                "parameters": subst_parameters,
            }
        )

    if arg.categories > 1 or arg.invariant:
        sitemodel_parameters = []
        if arg.categories > 1:
            sitemodel_parameters.append("&sitemodel.alpha")
        if arg.invariant:
            sitemodel_parameters.append("$sitemodel.proportions")
        list_opt.append(
            {
                "id": "opt.sitemodel",
                "type": "optimizer",
                "algorithm": "brent",
                "model": f"&{joint}",
                "parameters": sitemodel_parameters,
            }
        )

    meta["list"] = list_opt
    return meta


def create_tree_logger(id_, tree_id, arg):
    return {
        "id": id_,
        "type": "logger",
        "models": f"&{tree_id}",
        "tree": True,
        "internal": False,
        "file": arg.stem + ".trees",
    }


def build_optimizer(arg):
    json_list = []
    taxa = create_taxa("taxa", arg)

    if arg.map:
        joint_dic = create_evolution_joint(taxa, arg)
    else:
        joint_dic = create_tree_likelihood("treelikelihood", taxa, "alignment", arg)

    json_list.append(joint_dic)

    json_dic = {dic["id"]: dic for dic in json_list}

    json_dic["physher"] = []
    optimizer = create_meta_optimizer(joint_dic["id"], arg)
    json_dic["physher"].append(optimizer)

    tree_logger = create_tree_logger("treelogger", "tree", arg)
    json_dic["physher"].append(tree_logger)

    if arg.model == "GTR":
        json_dic["physher"].append(
            {"id": "log1", "type": "logger", "models": ["$frequencies", "$rates"]}
        )
    if arg.model == "HKY":
        json_dic["physher"].append(
            {"id": "log.kappa", "type": "logger", "parameters": ["&kappa"]}
        )
        json_dic["physher"].append(
            {"id": "log.frequencies", "type": "logger", "models": ["$frequencies"]}
        )

    if arg.categories > 1:
        json_dic["physher"].append(
            {"id": "log.alpha", "type": "logger", "parameters": ["&sitemodel.alpha"]}
        )
    if arg.invariant:
        json_dic["physher"].append(
            {"id": "log.pinv", "type": "logger", "models": ["$sitemodel.proportions"]}
        )
    return json_dic
