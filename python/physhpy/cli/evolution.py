import argparse
import logging
import os
import re
import sys

from dendropy import Tree

from .utils import (
    convert_date_to_real,
    extract_taxa,
    read_dates_from_csv,
    read_fasta_sequences,
)

logger = logging.getLogger(__name__)


def create_evolution_parser(parser):
    parser.add_argument("-i", "--input", required=False, help="""alignment file""")
    parser.add_argument("-t", "--tree", required=True, help="""tree file""")
    parser.add_argument(
        "-m",
        "--model",
        choices=["JC69", "HKY", "GTR"],
        default="JC69",
        help="""substitution model [default: %(default)s]""",
    )
    parser.add_argument(
        "-I",
        "--invariant",
        action="store_true",
        help="""include a proportion of invariant sites""",
    )
    parser.add_argument(
        "-f",
        "--frequencies",
        help="""frequencies""",
    )
    parser.add_argument(
        "-C",
        "--categories",
        metavar="C",
        type=int,
        default=1,
        help="""number of rate categories [default: %(default)s]""",
    )
    parser.add_argument(
        "--site_distr",
        choices=["gamma", "weibull"],
        default="weibull",
        help="""distribution for site model""",
    )
    parser.add_argument(
        "--brlenspr",
        choices=["exponential", "gammadir"],
        default="exponential",
        help="""prior on branch lengths of an unrooted tree [default: %(default)s]""",
    )
    parser.add_argument(
        "--clock",
        choices=["strict", "ucln"],
        help="""type of clock""",
    )
    parser.add_argument(
        "--clockpr",
        default="ctmcscale",
        type=lambda x: distribution_type(x, ("exponential", "ctmcscale")),
        help="""prior on substitution rate [default: %(default)s]""",
    )
    parser.add_argument(
        "--heights_init",
        choices=["tree"],
        help="""initialize node heights using input tree file or
         root to tip regression""",
    )
    parser.add_argument("--rate", type=float, help="""fixed substitution rate""")
    parser.add_argument(
        "--rate_init",
        type=lambda x: str_or_float(x, "regression"),
        help="""initialize substitution rate using 'regression' or with a value""",
    )
    parser.add_argument(
        "--dates",
        type=zero_or_path,
        help="""a csv file or 0 for contemporaneous taxa""",
    )
    parser.add_argument(
        "--date_format",
        default=None,
        help="""format of the date (yyyy/MM/dd or dd/MM/yyyy or dd-MM-yyyy)""",
    )
    parser.add_argument(
        "--date_regex",
        default=None,
        help="""regular expression to capture sampling date in sequence names""",
    )
    parser.add_argument(
        "--genetic_code",
        type=int,
        help="""index of genetic code""",
    )
    parser.add_argument(
        "--keep", action="store_true", help="""use branch length as starting values"""
    )
    parser.add_argument(
        "--use_path",
        action="store_true",
        help="""specify the alignment path instead of embedding it in the JSON file""",
    )
    parser.add_argument(
        "--use_ambiguities",
        action="store_true",
        help="""use nucleotide ambiguity codes""",
    )
    parser.add_argument(
        "--use_tip_states",
        action="store_true",
        help="""use tip states instead of tip partials""",
    )
    parser.add_argument(
        "--include_jacobian",
        action="store_true",
        help="""include Jacobian of the node height transform""",
    )

    parser = add_coalescent(parser)

    parser.add_argument(
        "--grid",
        type=int,
        help="""number of grid points (number of segments) for skygrid and BDSK""",
    )

    return parser


def add_coalescent(parser):
    parser.add_argument(
        "--coalescent",
        choices=["constant", "exponential", "skyride", "skygrid"],
        default=None,
        help="""type of coalescent""",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        help="""a cutoff for skygrid""",
    )
    parser.add_argument(
        "--time-aware",
        action="store_true",
        help="""use time aware GMRF""",
    )
    return parser


def zero_or_path(arg):
    if arg == "0":
        return 0
    elif arg is not None and not os.path.exists(arg):
        raise argparse.ArgumentTypeError(
            "invalid choice (choose from 0 or a path to a text file)"
        )
    else:
        return arg


def str_or_float(arg, choices):
    """Used by argparse when the argument can be either a number or a string
    from a prespecified list of options."""
    try:
        return float(arg)
    except ValueError:
        if (isinstance(choices, tuple) and arg in choices) or choices == arg:
            return arg
        else:
            if isinstance(choices, str):
                choices = (choices,)
            message = "'" + "','".join(choices) + '"'
            raise argparse.ArgumentTypeError(
                "invalid choice (choose from a number or " + message + ")"
            )


def distribution_type(arg, choices):
    """Used by argparse for specifying distributions with optional
    parameters."""
    res = arg.split("(")
    if (isinstance(choices, tuple) and res[0] in choices) or res[0] == choices:
        return arg
    else:
        if isinstance(choices, tuple):
            message = "'" + "','".join(choices) + '"'
        else:
            message = "'" + choices + "'"
        raise argparse.ArgumentTypeError(
            "invalid choice (choose from a number or " + message + ")"
        )


def create_tree_model(id_: str, taxa: dict, arg):
    tree_model = {
        "id": id_,
        "type": "tree",
    }
    if arg.clock is not None:
        dates = [taxon["attributes"]["date"] for taxon in taxa]
        offset = max(dates) - min(dates)

        tree_model["heights"] = "tree.heights"
        # tree_model["reparam"] = "tree.scalers"
        tree_model["time"] = True
        tree_model["dates"] = {
            taxon["id"]: taxon["attributes"]["date"] for taxon in taxa
        }
        tree_model["ratios"] = {
            "id": f"{id_}.ratios",
            "type": "parameter",
            "dimension": len(taxa) - 2,
            "values": [0.1],
            "lower": 0,
            "upper": 1,
        }
        tree_model["root_height"] = {
            "id": f"{id_}.root_height",
            "type": "parameter",
            "value": offset + 1.0,
            "lower": offset,
        }
    else:
        tree_model["branch_lengths"] = {
            "id": f"{id_}.distances",
            "type": "parameter",
            "dimension": 2 * len(taxa) - 3,
            "values": [0.1],
            "lower": 0
            # "values": [0.1]*(2*len(taxa)-3)
        }

    tree_format = "newick"
    with open(arg.tree, "r") as fp:
        if next(fp).upper().startswith("#NEXUS"):
            tree_format = "nexus"
    if tree_format == "nexus":
        tree = Tree.get(
            path=arg.tree,
            schema=tree_format,
            tree_offset=0,
            preserve_underscores=True,
        )
        newick = str(tree) + ";"
    else:
        with open(arg.tree, "r") as fp:
            newick = fp.read()
            newick = newick.strip()

    tree_model["newick"] = newick

    return tree_model


def create_poisson_tree_likelihood(id_, taxa, arg):
    tree_id = "tree"
    tree_model = create_tree_model(tree_id, taxa, arg)
    branch_model = create_branch_model("branchmodel", tree_id, len(taxa["taxa"]), arg)

    treelikelihood_model = {
        "id": id_,
        "type": "PoissonTreeLikelihood",
        "tree_model": tree_model,
        "branch_model": branch_model,
    }

    return treelikelihood_model


def create_tree_likelihood_single(
    id_, tree_model, branch_model, substitution_model, site_model, site_pattern
):

    treelikelihood_model = {
        "id": id_,
        "type": "treelikelihood",
        "tree": tree_model,
        "sitemodel": site_model,
        "substitutionmodel": substitution_model,
        "sitepattern": site_pattern,
    }
    if branch_model is not None:
        treelikelihood_model["branchmodel"] = branch_model

    return treelikelihood_model


def create_tree_likelihood(id_, taxa, alignment, arg):
    site_pattern = create_site_pattern("patterns", arg)
    site_model = create_site_model("sitemodel", arg)
    substitution_model = create_substitution_model("substmodel", arg.model, arg)
    tree_id = "tree"
    tree_model = create_tree_model(tree_id, taxa, arg)

    treelikelihood_model = {
        "id": id_,
        "type": "treelikelihood",
        "tree": tree_model,
        "sitemodel": site_model,
        "substitutionmodel": substitution_model,
        "sitepattern": site_pattern,
    }

    if arg.use_tip_states:
        treelikelihood_model["tipstates"] = True

    if arg.include_jacobian:
        treelikelihood_model["include_jacobian"] = True

    if arg.clock is not None:
        treelikelihood_model["branchmodel"] = create_branch_model(
            "branchmodel", tree_id, len(taxa), arg, arg.rate_init
        )

    return treelikelihood_model


def create_site_model(id_, arg, w=None):
    site_model = {"id": id_, "type": "sitemodel"}
    distribution = None
    if arg.categories > 1:
        distribution = {
            "distribution": arg.site_distr,
            "categories": arg.categories,
            "parameters": {
                "shape": {
                    "id": f"{id_}.alpha",
                    "type": "parameter",
                    "value": 1,
                    "lower": 0,
                }
            },
        }
        if arg.invariant:
            distribution["proportions"] = {
                "id": f"{id_}.proportions",
                "type": "simplex",
                "dimension": 2,
            }
    elif arg.invariant:
        distribution = {
            "type": "distribution",
            "distribution": "discrete",
            "categories": 2,
            "proportions": {
                "id": f"{id_}.proportions",
                "type": "simplex",
                "dimension": 2,
                "values": [0.5, 0.5],
            },
        }
    if distribution is not None:
        site_model["distribution"] = distribution

    if arg.model == "SRD06":
        site_model["mu"] = {
            "id": f"{id_}.mu",
            "type": "parameter",
            "value": 1,
            "lower": 0,
        }
    return site_model


def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False
    except TypeError:
        return False


def create_branch_model(id_, tree_id, taxa_count, arg, rate_init=None):
    branch_model = {"id": id_, "type": "branchmodel", "tree": "&tree"}
    if arg.clock == "strict":
        branch_model["model"] = "strict"
        branch_model["rate"] = {
            "id": "rate",
            "type": "parameter",
            "value": 0.003,
            "lower": 0.0,
        }

    if rate_init is not None:
        branch_model["rate"] = rate_init
    return branch_model


def create_substitution_model(id_, model, arg):
    subst_model = {"id": id_, "type": "substitutionmodel"}
    subst_model["model"] = model.lower()
    subst_model["frequencies"] = {
        "id": "frequencies",
        "type": "simplex",
        "values": [0.25] * 4,
    }
    if model == "JC69":
        subst_model["frequencies"]["!PHYSHER_FIXED"] = True

    if model in ("JC69", "HKY", "GTR"):
        subst_model["datatype"] = "nucleotide"

    if model == "GTR":
        subst_model["rates"] = {"id": "rates", "type": "simplex", "values": [1 / 6] * 6}
    elif model == "HKY":
        subst_model["rates"] = {
            "kappa": {"id": "kappa", "type": "parameter", "value": 1, "lower": 0}
        }
    return subst_model


def create_site_pattern(id_, arg):
    site_pattern = {
        "id": id_,
        "type": "sitepattern",
        "datatype": "nucleotide",
        "alignment": {
            "id": "seqs",
            "type": "alignment",
        },
    }
    if arg.use_path:
        site_pattern["alignment"]["file"] = arg.input
    else:
        sequences = read_fasta_sequences(arg.input)
        site_pattern["alignment"]["sequences"] = {
            seq.taxon: seq.sequence for seq in sequences
        }
    return site_pattern


def create_taxa(id_, arg):
    if arg.input is not None:
        alignment = read_fasta_sequences(arg.input)
        taxa_list = [{"id": sequence.taxon} for sequence in alignment]
    else:
        taxa_list = extract_taxa(arg.tree)

    if arg.clock is not None:
        if arg.dates == 0:
            for taxon in taxa_list:
                taxon["attributes"] = {"date": 0.0}
        elif arg.dates is not None:
            dates = read_dates_from_csv(arg.dates, arg.date_format)
            for taxon in taxa_list:
                taxon["attributes"] = {"date": dates[taxon["id"]]}
        else:
            regex_date = r"_(\d+\.?\d*)$"
            if arg.date_regex is not None:
                regex_date = arg.date_regex
            regex = re.compile(regex_date)

            if arg.date_format is not None:
                res = re.split(r"[/-]", arg.date_format)
                yy = res.index("yyyy") + 1
                MM = res.index("MM") + 1
                dd = res.index("dd") + 1

            for taxon in taxa_list:
                res = re.search(regex, taxon["id"])
                if res is None:
                    logger.error(
                        f" Could not extract date from {taxon['id']}"
                        f" - regular expression used: {regex_date}"
                    )
                    sys.exit(1)
                if len(res.groups()) > 1:
                    taxon["attributes"] = {
                        "date": convert_date_to_real(
                            int(res[dd]), int(res[MM]), int(res[yy])
                        )
                    }
                else:
                    taxon["attributes"] = {"date": float(res.group(1))}
    return taxa_list


def create_coalesent(id_, tree_id, theta_id, arg, **kwargs):
    if arg.coalescent == "constant":
        coalescent = {
            "id": id_,
            "type": "ConstantCoalescentModel",
            "theta": theta_id,
            "tree_model": tree_id,
        }
    elif arg.coalescent == "exponential":
        coalescent = {
            "id": id_,
            "type": "ExponentialCoalescentModel",
            "theta": theta_id,
            "growth": kwargs.get("growth"),
            "tree_model": tree_id,
        }
    elif arg.coalescent == "skygrid":
        coalescent = {
            "id": id_,
            "type": "PiecewiseConstantCoalescentGridModel",
            "theta": theta_id,
            "tree_model": tree_id,
            "cutoff": arg.cutoff,
        }
    elif arg.coalescent == "skyride":
        coalescent = {
            "id": id_,
            "type": "PiecewiseConstantCoalescentModel",
            "theta": theta_id,
            "tree_model": tree_id,
        }

    return coalescent


def create_substitution_model_priors(substmodel_id, model):
    joint_list = []
    if model == "HKY" or model == "GTR":
        joint_list.append(
            {
                "id": "priorfreqs",
                "type": "distribution",
                "distribution": "dirichlet",
                "x": "&frequencies",
            }
        )

        if model == "GTR":
            joint_list.append(
                {
                    "id": "priorrates",
                    "type": "distribution",
                    "distribution": "dirichlet",
                    "x": "&rates",
                }
            )
        else:
            joint_list.append(
                {
                    "id": "priorrates",
                    "type": "distribution",
                    "distribution": "lognormal",
                    "parameters": {
                        "mu": {
                            "id": "mean.kappa",
                            "type": "parameter",
                            "value": 1,
                            "!PHYSHER_FIXED": True,
                        },
                        "sigma": {
                            "id": "sigma.kappa",
                            "type": "parameter",
                            "value": 1.25,
                            "lower": 0.0,
                            "!PHYSHER_FIXED": True,
                        },
                    },
                    "x": ["&kappa"],
                }
            )
    return joint_list


def parse_distribution(desc):
    res = desc.split("(")
    if len(res) == 1:
        return res, None
    else:
        return res[0], list(map(float, res[1][:-1].split(",")))


def create_clock_prior(arg):
    tree_id = "tree"
    prior_list = []
    if arg.clock == "strict" and arg.rate is None:
        if arg.clockpr == "ctmcscale":
            prior_list.append(
                {
                    "id": "priorrate",
                    "type": "distribution",
                    "distribution": "ctmcscale",
                    "x": "&rate",
                    "tree": f"&{tree_id}",
                }
            )
        elif arg.clockpr.startswith("exponential"):
            name, params = parse_distribution(arg.clockpr)
            if params is None:
                rate = 1000.0
            else:
                rate = params[0]
            prior_list.append(
                {
                    "id": "prior.clock",
                    "type": "distribution",
                    "distribution": "exponential",
                    "parameters": {
                        "mean": {
                            "id": "prior.clock.mean",
                            "type": "parameter",
                            "value": rate,
                            "lower": 0,
                        }
                    },
                    "x": "&rate",
                },
            )
    return prior_list


def create_evolution_priors(taxa, arg):
    tree_id = "tree"
    joint_list = []

    if arg.coalescent is not None:
        joint_list.extend(create_time_tree_prior(taxa, arg))

    if arg.clock is not None:
        joint_list.extend(create_clock_prior(arg))
    else:
        if arg.brlenspr == "exponential":
            joint_list.append(
                {
                    "id": "prior.distances",
                    "type": "distribution",
                    "distribution": "exponential",
                    "parameters": {
                        "mean": {
                            "id": "prior.distances.mean",
                            "type": "parameter",
                            "value": 0.1,
                            "lower": 0,
                            "!PHYSHER_FIXED": True,
                        }
                    },
                    "x": f"%{tree_id}.distances",
                },
            )

    joint_list.extend(create_substitution_model_priors("substmodel", arg.model))

    if arg.categories > 1:
        sitemodel_id = "sitemodel"
        joint_list.append(
            {
                "id": "prior.sitemodel",
                "type": "distribution",
                "distribution": "exponential",
                "parameters": {
                    "mean": {
                        "id": "prior.sitemodel.mean",
                        "type": "parameter",
                        "value": 1.0,
                        "lower": 0,
                    }
                },
                "x": f"&{sitemodel_id}.alpha",
            },
        )
    return joint_list


def create_time_tree_prior(taxa, arg):
    tree_id = "tree"
    joint_list = []
    if arg.coalescent is not None:
        coalescent = {
            "id": "coalescent",
            "type": "coalescent",
            "model": arg.coalescent,
            "parameters": {},
            "tree": f"&{tree_id}",
        }
        coalescent_id = "coalescent"
        if arg.coalescent in ("constant", "exponential"):
            coalescent["parameters"]["n0"] = {
                "id": "n0",
                "type": "parameter",
                "value": 3,
                "lower": 0.00,
                "upper": "infinity",
            }

            joint_list.append(
                {
                    "id": f"{coalescent_id}.theta.prior",
                    "type": "distribution",
                    "distribution": "oneonx",
                    "x": "&n0",
                }
            )
        if arg.coalescent == "exponential":
            coalescent["parameters"]["growth"] = {
                "id": "n0",
                "type": "parameter",
                "value": 1,
            }
        elif arg.coalescent == "skygrid":
            coalescent["parameters"]["theta"] = {
                "id": "thetas",
                "type": "parameter",
                "values": [3.0] * arg.grid,
            }
            coalescent["parameters"]["cutoff"] = (arg.cutoff,)
        elif arg.coalescent == "skyride":
            coalescent["parameters"]["theta"] = {
                "id": "thetas",
                "type": "parameter",
                "values": [3.0] * (len(taxa["taxa"]) - 1),
            }

    return [coalescent] + joint_list


def create_evolution_joint(taxa, arg):
    joint_list = [
        create_tree_likelihood("treelikelihood", taxa, arg),
    ] + create_evolution_priors(taxa, arg)

    joint_dic = {
        "id": "joint",
        "type": "compound",
        "distributions": joint_list,
    }
    return joint_dic
