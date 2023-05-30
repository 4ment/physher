// Copyright 2016-2022 Mathieu Fourment.
// physher is free software under the GPLv2; see LICENSE file for details.

#include "phycpp/physher.hpp"

#include <algorithm>
#include <iostream>

extern "C" {
#include "phyc/ctmcscale.h"
#include "phyc/distmodel.h"
#include "phyc/treetransform.h"
}

void TreeModelInterface::InitializeMap(const std::vector<std::string> &taxa) {
    nodeMap_.resize(nodeCount_);
    for (size_t i = 0; i < nodeCount_; i++) {
        Node *node = Tree_node(tree_, i);
        if (Node_isleaf(node)) {
            std::string str(Node_name(node));
            auto it = find(taxa.begin(), taxa.end(), str);
            if (it != taxa.end()) {
                nodeMap_[node->id] = it - taxa.begin();
            } else {
                std::cerr << "Could not find taxon " << str << " in taxon list"
                          << std::endl;
            }
        } else {
            nodeMap_[node->id] = node->class_id + tipCount_;
        }
    }
}

UnRootedTreeModelInterface::UnRootedTreeModelInterface(
    const std::string &newick, const std::vector<std::string> &taxa) {
    char **taxa_p = new char *[taxa.size()];
    for (size_t i = 0; i < taxa.size(); i++) {
        taxa_p[i] = const_cast<char *>(taxa[i].c_str());
    }
    model_ = new_TreeModel_from_newick(newick.c_str(), taxa_p, nullptr);
    delete[] taxa_p;

    tree_ = reinterpret_cast<Tree *>(model_->obj);
    nodeCount_ = Tree_node_count(tree_);
    tipCount_ = Tree_tip_count(tree_);

    InitializeMap(taxa);
    parameterCount_ = nodeCount_ - 2;
}

void UnRootedTreeModelInterface::SetParameters(const double *parameters) {
    Node **nodes = Tree_nodes(tree_);
    Node *root = Tree_root(tree_);
    for (size_t i = 0; i < nodeCount_; i++) {
        Node *node = nodes[i];
        if (node == root || node == root->left || node == root->right) continue;
        Node_set_distance(node, parameters[nodeMap_[node->id]]);
    }
    if (Node_isleaf(root->right)) {
        Node_set_distance(root->left, parameters[nodeMap_[root->right->id]]);
    } else {
        Node_set_distance(root->left, parameters[nodeMap_[root->left->id]]);
    }
}

ReparameterizedTimeTreeModelInterface::ReparameterizedTimeTreeModelInterface(
    const std::string &newick, const std::vector<std::string> &taxa,
    const std::vector<double> dates, TreeTransformFlags transform) {
    char **taxa_p = new char *[taxa.size()];
    for (size_t i = 0; i < taxa.size(); i++) {
        taxa_p[i] = const_cast<char *>(taxa[i].c_str());
    }
    model_ = new_TimeTreeModel_from_newick(newick.c_str(), taxa_p, dates.data());
    delete[] taxa_p;
    int int_transform =
        static_cast<std::underlying_type<TreeLikelihoodGradientFlags>::type>(transform);
    TreeModel_set_transform(model_, int_transform);

    tree_ = reinterpret_cast<Tree *>(model_->obj);
    nodeCount_ = Tree_node_count(tree_);
    tipCount_ = Tree_tip_count(tree_);

    InitializeMap(taxa);
    transformModel_ = reinterpret_cast<Model *>(model_->data);
    parameterCount_ = tipCount_ - 1;
}

void ReparameterizedTimeTreeModelInterface::SetParameters(const double *parameters) {
    TreeTransform *tt = reinterpret_cast<TreeTransform *>(transformModel_->obj);
    Parameters_set_values_quietly(tt->parameters, parameters);
    transformModel_->listeners->fire(transformModel_->listeners, transformModel_, -1);
}

void ReparameterizedTimeTreeModelInterface::GetNodeHeights(double *heights) {
    Tree_update_heights(tree_);
    Node **nodes = Tree_nodes(tree_);
    for (size_t i = 0; i < nodeCount_; i++) {
        if (!Node_isleaf(nodes[i])) {
            heights[nodes[i]->class_id] = Node_height(nodes[i]);
        }
    }
}

void ReparameterizedTimeTreeModelInterface::GradientTransformJVP(
    double *gradient, const double *height_gradient) {
    Tree_update_heights(tree_);
    Tree_node_transform_jvp(tree_, height_gradient, gradient);
}

void ReparameterizedTimeTreeModelInterface::GradientTransformJVP(
    double *gradient, const double *height_gradient, const double *heights) {
    Tree_node_transform_jvp_with_heights(tree_, heights, height_gradient, gradient);
}

void ReparameterizedTimeTreeModelInterface::GradientTransformJacobian(
    double *gradient) {
    memset(gradient, 0.0, sizeof(double) * (tipCount_ - 1));
    Tree_update_heights(tree_);
    Tree_node_transform_jacobian_gradient(tree_, gradient);
}

double ReparameterizedTimeTreeModelInterface::TransformJacobian() {
    Tree_update_heights(tree_);
    TreeTransform *tt = reinterpret_cast<TreeTransform *>(transformModel_->obj);
    return tt->log_jacobian(tt);
}

void BranchModelInterface::SetParameters(const double *parameters) {
    Parameters_set_values(branchModel_->rates, parameters);
}
void BranchModelInterface::GetParameters(double *parameters) {
    for (size_t i = 0; i < Parameters_count(branchModel_->rates); i++) {
        parameters[i] = Parameters_value(branchModel_->rates, i);
    }
}
void BranchModelInterface::SetRates(const double *rates) {
    Parameters_set_values(branchModel_->rates, rates);
}

StrictClockModelInterface::StrictClockModelInterface(double rate,
                                                     TreeModelInterface *treeModel) {
    Parameter *p = new_Parameter("clock.rate", rate, new_Constraint(0, INFINITY));
    p->model = MODEL_BRANCHMODEL;
    branchModel_ = new_StrictClock_with_parameter(treeModel->GetTree(), p);
    free_Parameter(p);
    model_ = new_BranchModel2("", branchModel_, treeModel->GetModel(), nullptr);
    parameterCount_ = 1;
}

void StrictClockModelInterface::SetRate(double rate) {
    Parameters_set_value(branchModel_->rates, 0, rate);
}

SimpleClockModelInterface::SimpleClockModelInterface(const std::vector<double> &rates,
                                                     TreeModelInterface *treeModel) {
    size_t nodeCount = treeModel->GetNodeCount();
    Parameters *ps = new_Parameters(nodeCount - 1);
    DiscreteParameter *map = new_DiscreteParameter(nodeCount);
    StringBuffer *buffer = new_StringBuffer(10);
    Tree *tree = treeModel->GetTree();
    for (int i = 0; i < nodeCount; i++) {
        Node *node = Tree_node(tree, i);
        if (!Node_isroot(node)) {
            size_t index = Node_isleaf(node)
                               ? node->class_id
                               : node->class_id + treeModel->GetTipCount();
            map->values[Node_id(node)] = index;
            StringBuffer_set_string(buffer, "clock_rates");
            StringBuffer_append_format(buffer, ".%d", index);
            Parameter *p =
                new_Parameter(buffer->c, rates[index], new_Constraint(0, INFINITY));
            p->id = Node_id(node);
            p->model = MODEL_BRANCHMODEL;
            Parameters_move(ps, p);
        }
    }
    parameterCount_ = nodeCount - 1;
    free_StringBuffer(buffer);
    branchModel_ = new_DiscreteClock_with_parameters(treeModel->GetTree(), ps, map);
    free_Parameters(ps);
    model_ =
        new_BranchModel2("simple_clock", branchModel_, treeModel->GetModel(), nullptr);
}

JC69Interface::JC69Interface() {
    Simplex *frequencies_simplex = new_Simplex("jc69_frequency_simplex", 4);
    Model *frequencies_model =
        new_SimplexModel("jc69_frequencies", frequencies_simplex);
    dataType_ = new NucleotideDataTypeInterface();
    substModel_ = new_JC69(frequencies_simplex);
    free_DataType(substModel_->datatype);
    substModel_->datatype = dataType_->dataType_;
    substModel_->datatype->ref_count++;
    model_ = new_SubstitutionModel2("jc69", substModel_, frequencies_model, nullptr);
    frequencies_model->free(frequencies_model);
    parameterCount_ = 0;
}

HKYInterface::HKYInterface(double kappa, const std::vector<double> &frequencies) {
    Parameter *kappa_parameter =
        new_Parameter("hky_kappa", kappa, new_Constraint(0, INFINITY));
    Simplex *frequencies_simplex = new_Simplex_with_values(
        "hky_frequency_simplex", frequencies.data(), frequencies.size());
    Model *frequencies_model = new_SimplexModel("hky_frequencies", frequencies_simplex);
    dataType_ = new NucleotideDataTypeInterface();
    substModel_ = new_HKY_with_parameters(frequencies_simplex, kappa_parameter);
    free_DataType(substModel_->datatype);
    substModel_->datatype = dataType_->dataType_;
    substModel_->datatype->ref_count++;
    model_ = new_SubstitutionModel2("hky", substModel_, frequencies_model, nullptr);
    free_Parameter(kappa_parameter);
    frequencies_model->free(frequencies_model);
    parameterCount_ = 5;
}

void HKYInterface::SetKappa(double kappa) {
    Parameters_set_value(substModel_->rates, 0, kappa);
}

void HKYInterface::SetFrequencies(const double *frequencies) {
    substModel_->simplex->set_values(substModel_->simplex, frequencies);
}

void HKYInterface::SetParameters(const double *parameters) {
    SetKappa(parameters[1]);
    substModel_->simplex->set_values(substModel_->simplex, parameters + 1);
}

GTRInterface::GTRInterface(const std::vector<double> &rates,
                           const std::vector<double> &frequencies) {
    Parameters *rates_parameters = nullptr;
    Model *rates_model = nullptr;
    // simplex
    if (rates.size() == 6) {
        Simplex *rates_simplex =
            new_Simplex_with_values("gtr_rate_simplex", rates.data(), rates.size());
        rates_model = new_SimplexModel("gtr_rates", rates_simplex);
        parameterCount_ = 10;
    } else {
        rates_parameters = new_Parameters(5);
        for (auto rate : rates) {
            Parameters_move(
                rates_parameters,
                new_Parameter("gtr_rates", rate, new_Constraint(0, INFINITY)));
        }
        parameterCount_ = 9;
    }
    Simplex *frequencies_simplex = new_Simplex_with_values(
        "gtr_frequency_simplex", frequencies.data(), frequencies.size());
    Model *frequencies_model = new_SimplexModel("gtr_frequencies", frequencies_simplex);

    dataType_ = new NucleotideDataTypeInterface();
    const char *assigments[5] = {"AC", "AG", "AT", "CG", "CT"};
    Simplex *rates_simplex = rates_model == nullptr
                                 ? nullptr
                                 : reinterpret_cast<Simplex *>(rates_model->obj);
    substModel_ =
        SubstitutionModel_factory("gtr", dataType_->dataType_,
                                  reinterpret_cast<Simplex *>(frequencies_model->obj),
                                  rates_simplex, rates_parameters, assigments);
    free_DataType(substModel_->datatype);
    substModel_->datatype = dataType_->dataType_;
    substModel_->datatype->ref_count++;
    model_ = new_SubstitutionModel2("gtr", substModel_, frequencies_model, rates_model);
    free_Parameters(rates_parameters);
    frequencies_model->free(frequencies_model);
    if (rates_model != nullptr) {
        rates_model->free(rates_model);
    }
}

void GTRInterface::SetRates(const double *rates) {
    if (substModel_->rates == nullptr) {
        substModel_->rates_simplex->set_values(substModel_->rates_simplex, rates);
    } else {
        Parameters_set_values(substModel_->rates, rates);
    }
}
void GTRInterface::SetFrequencies(const double *frequencies) {
    substModel_->simplex->set_values(substModel_->simplex, frequencies);
}
void GTRInterface::SetParameters(const double *parameters) {
    Parameters_set_values(substModel_->rates, parameters);
    substModel_->simplex->set_values(substModel_->simplex, parameters + 3);
}

GeneralSubstitutionModelInterface::GeneralSubstitutionModelInterface(
    DataTypeInterface *dataType, const std::vector<double> &rates,
    const std::vector<double> &frequencies, const std::vector<unsigned> &structure,
    bool normalize) {
    Parameters *rate_parameters = new_Parameters(rates.size());
    for (auto rate : rates) {
        Parameters_move(rate_parameters, new_Parameter("subst_rates", rate,
                                                       new_Constraint(0, INFINITY)));
    }

    parameterCount_ = rates.size();
    Simplex *frequencies_simplex = new_Simplex_with_values(
        "subst_frequency_simplex", frequencies.data(), frequencies.size());
    Model *frequencies_model =
        new_SimplexModel("subst_frequencies", frequencies_simplex);

    DiscreteParameter *dp =
        new_DiscreteParameter_with_values(structure.data(), structure.size());
    Model *mdp = new_DiscreteParameterModel("structure", dp);

    substModel_ = new_GeneralModel_with_parameters(
        dataType->dataType_, reinterpret_cast<DiscreteParameter *>(mdp->obj),
        rate_parameters, frequencies_simplex, -1, normalize);
    model_ = new_SubstitutionModel3("substmodel", substModel_, frequencies_model,
                                    nullptr, mdp);

    free_Parameters(rate_parameters);
    frequencies_model->free(frequencies_model);
    dataType_ = dataType;
}

void GeneralSubstitutionModelInterface::SetRates(const double *rates) {
    if (substModel_->rates == nullptr) {
        substModel_->rates_simplex->set_values(substModel_->rates_simplex, rates);
    } else {
        Parameters_set_values(substModel_->rates, rates);
    }
}

void GeneralSubstitutionModelInterface::SetFrequencies(const double *frequencies) {
    substModel_->simplex->set_values(substModel_->simplex, frequencies);
}

void GeneralSubstitutionModelInterface::SetParameters(const double *parameters) {
    Parameters_set_values(substModel_->rates, parameters);
    substModel_->simplex->set_values(substModel_->simplex, parameters + 3);
}

void SiteModelInterface::SetMu(double mu) { Parameter_set_value(siteModel_->mu, mu); }

ConstantSiteModelInterface::ConstantSiteModelInterface(std::optional<double> mu) {
    siteModel_ = new_SiteModel_with_parameters(
        nullptr, nullptr, 1, DISTRIBUTION_UNIFORM, false, QUADRATURE_QUANTILE_MEDIAN);
    if (mu.has_value()) {
        Parameter *mu_param = new_Parameter("mu", *mu, new_Constraint(0, INFINITY));
        SiteModel_set_mu(siteModel_, mu_param);
        free_Parameter(mu_param);
    }
    model_ = new_SiteModel2("sitemodel", siteModel_, nullptr);
    parameterCount_ = (mu.has_value() ? 1 : 0);
}

void ConstantSiteModelInterface::SetParameters(const double *parameters) {
    if (siteModel_->mu != nullptr) {
        Parameter_set_value(siteModel_->mu, parameters[0]);
    }
}
void ConstantSiteModelInterface::GetParameters(double *parameters) {
    if (siteModel_->mu != nullptr) {
        parameters[0] = Parameter_value(siteModel_->mu);
    }
}

InvariantSiteModelInterface::InvariantSiteModelInterface(double proportionInvariant,
                                                         std::optional<double> mu) {
    Simplex *propInvSimplex = new_Simplex("pinv", 2);
    Model *propInvModel = new_SimplexModel("pinv.model", propInvSimplex);
    Parameters_at(propInvSimplex->parameters, 0)->model = MODEL_SITEMODEL;
    Parameters_set_value(propInvSimplex->parameters, 0, proportionInvariant);

    siteModel_ =
        new_SiteModel_with_parameters(nullptr, propInvSimplex, 1, DISTRIBUTION_DISCRETE,
                                      true, QUADRATURE_QUANTILE_MEDIAN);
    if (mu.has_value()) {
        Parameter *mu_param = new_Parameter("mu", *mu, new_Constraint(0, INFINITY));
        SiteModel_set_mu(siteModel_, mu_param);
        free_Parameter(mu_param);
    }

    model_ = new_SiteModel2("sitemodel", siteModel_, propInvModel);
    parameterCount_ = 1 + siteModel_->mu != nullptr;
}

void InvariantSiteModelInterface::SetProportionInvariant(double value) {
    Parameters_set_value(siteModel_->proportions->parameters, 0, value);
}

void InvariantSiteModelInterface::SetParameters(const double *parameters) {
    size_t shift = 0;
    if (siteModel_->proportions != nullptr) {
        Parameters_set_value(siteModel_->proportions->parameters, 0, parameters[shift]);
        shift++;
    }

    if (siteModel_->mu != nullptr) {
        Parameter_set_value(siteModel_->mu, parameters[shift]);
    }
}

void InvariantSiteModelInterface::GetParameters(double *parameters) {
    size_t shift = 0;
    if (siteModel_->proportions != nullptr) {
        parameters[shift] = Parameters_value(siteModel_->proportions->parameters, 0);
        shift++;
    }
    if (siteModel_->mu != nullptr) {
        parameters[shift] = Parameter_value(siteModel_->mu);
    }
}

DiscretizedSiteModelInterface::DiscretizedSiteModelInterface(
    distribution_t distribution, double shape, size_t categories,
    std::optional<double> proportionInvariant, std::optional<double> mu) {
    Parameter *shape_parameter =
        new_Parameter("sitemodel.parameter", shape, new_Constraint(0, INFINITY));
    Parameters *params = new_Parameters(1);
    Parameters_move(params, shape_parameter);

    Simplex *propInvSimplex = nullptr;
    Model *propInvModel = nullptr;
    if (proportionInvariant.has_value()) {
        propInvSimplex = new_Simplex("pinv", 2);
        propInvModel = new_SimplexModel("pinv.model", propInvSimplex);
        Parameters_at(propInvSimplex->parameters, 0)->model = MODEL_SITEMODEL;
        Parameters_set_value(propInvSimplex->parameters, 0, *proportionInvariant);
    }

    siteModel_ = new_SiteModel_with_parameters(
        params, propInvSimplex, categories, distribution,
        proportionInvariant.has_value(), QUADRATURE_QUANTILE_MEDIAN);
    if (mu.has_value()) {
        Parameter *mu_param = new_Parameter("mu", *mu, new_Constraint(0, INFINITY));
        SiteModel_set_mu(siteModel_, mu_param);
        free_Parameter(mu_param);
    }
    siteModel_->epsilon = 1.e-6;

    model_ = new_SiteModel2("sitemodel", siteModel_, propInvModel);
    free_Parameters(params);
    parameterCount_ = Parameters_count(siteModel_->rates) + siteModel_->mu != nullptr;
    if (propInvModel != nullptr) {
        propInvModel->free(propInvModel);
        parameterCount_ += propInvSimplex->K - 1;
    }
}

void DiscretizedSiteModelInterface::SetParameter(double parameter) {
    Parameters_set_value(siteModel_->rates, 0, parameter);
}

void DiscretizedSiteModelInterface::SetProportionInvariant(double value) {
    Parameters_set_value(siteModel_->proportions->parameters, 0, value);
}

void DiscretizedSiteModelInterface::SetParameters(const double *parameters) {
    size_t shift = 0;
    if (Parameters_count(siteModel_->rates) > 0) {
        Parameters_set_value(siteModel_->rates, 0, parameters[0]);
        shift++;
    }
    if (siteModel_->proportions != nullptr) {
        Parameters_set_value(siteModel_->proportions->parameters, 0, parameters[shift]);
        shift++;
    }

    if (siteModel_->mu != nullptr) {
        Parameter_set_value(siteModel_->mu, parameters[shift]);
        shift++;
    }
}

void DiscretizedSiteModelInterface::GetParameters(double *parameters) {
    size_t shift = 0;
    if (Parameters_count(siteModel_->rates) > 0) {
        parameters[0] = Parameters_value(siteModel_->rates, 0);
        shift++;
    }
    if (siteModel_->proportions != nullptr) {
        parameters[shift] = Parameters_value(siteModel_->proportions->parameters, 0);
        shift++;
    }
    if (siteModel_->mu != nullptr) {
        parameters[shift] = Parameter_value(siteModel_->mu);
        shift++;
    }
}

void WeibullSiteModelInterface::SetShape(double shape) {
    Parameters_set_value(siteModel_->rates, 0, shape);
}

void GammaSiteModelInterface::SetShape(double shape) {
    Parameters_set_value(siteModel_->rates, 0, shape);
}

void GammaSiteModelInterface::SetEpsilon(double epsilon) {
    siteModel_->epsilon = epsilon;
}

TreeLikelihoodInterface::TreeLikelihoodInterface(
    const std::vector<std::pair<std::string, std::string>> &alignment,
    TreeModelInterface *treeModel, SubstitutionModelInterface *substitutionModel,
    SiteModelInterface *siteModel, std::optional<BranchModelInterface *> branchModel,
    bool use_ambiguities, bool use_tip_states, bool include_jacobian)
    : treeModel_(treeModel),
      substitutionModel_(substitutionModel),
      siteModel_(siteModel) {
    branchModel_ = branchModel.has_value() ? *branchModel : nullptr;
    Sequences *sequences = new_Sequences(alignment.size());
    for (const auto &sequence : alignment) {
        Sequences_add(sequences,
                      new_Sequence(sequence.first.c_str(), sequence.second.c_str()));
    }
    sequences->datatype = substitutionModel->GetDataType()->dataType_;
    sequences->datatype->ref_count++;
    SitePattern *sitePattern = new_SitePattern(sequences);
    free_Sequences(sequences);
    Model *mbm = branchModel.has_value() ? branchModel_->GetModel() : nullptr;
    BranchModel *bm =
        branchModel.has_value() ? reinterpret_cast<BranchModel *>(mbm->obj) : nullptr;
    SingleTreeLikelihood *tlk = new_SingleTreeLikelihood(
        reinterpret_cast<Tree *>(treeModel_->GetManagedObject()),
        reinterpret_cast<SubstitutionModel *>(substitutionModel_->GetManagedObject()),
        reinterpret_cast<SiteModel *>(siteModel_->GetManagedObject()), sitePattern, bm,
        use_tip_states);
    model_ = new_TreeLikelihoodModel("treelike", tlk, treeModel_->GetModel(),
                                     substitutionModel_->GetModel(),
                                     siteModel_->GetModel(), mbm);

    tlk->include_jacobian = include_jacobian;
    RequestGradient();
}

TreeLikelihoodInterface::TreeLikelihoodInterface(
    const std::vector<std::string> &taxa, const std::vector<std::string> &attributes,
    TreeModelInterface *treeModel, SubstitutionModelInterface *substitutionModel,
    SiteModelInterface *siteModel, std::optional<BranchModelInterface *> branchModel,
    bool use_ambiguities, bool use_tip_states, bool include_jacobian)
    : treeModel_(treeModel),
      substitutionModel_(substitutionModel),
      siteModel_(siteModel) {
    branchModel_ = branchModel.has_value() ? *branchModel : nullptr;
    const char **taxaPtr = new const char *[taxa.size()];
    const char **attributesPtr = new const char *[taxa.size()];
    for (size_t i = 0; i < taxa.size(); i++) {
        taxaPtr[i] = const_cast<char *>(taxa[i].c_str());
        attributesPtr[i] = const_cast<char *>(attributes[i].c_str());
    }
    DataType *dataType = substitutionModel->GetDataType()->dataType_;
    SitePattern *sitePattern =
        new_AttributePattern(dataType, taxaPtr, attributesPtr, taxa.size());
    delete[] taxaPtr;
    delete[] attributesPtr;

    Model *mbm = branchModel.has_value() ? branchModel_->GetModel() : nullptr;
    BranchModel *bm =
        branchModel.has_value() ? reinterpret_cast<BranchModel *>(mbm->obj) : nullptr;
    SingleTreeLikelihood *tlk = new_SingleTreeLikelihood(
        reinterpret_cast<Tree *>(treeModel_->GetManagedObject()),
        reinterpret_cast<SubstitutionModel *>(substitutionModel_->GetManagedObject()),
        reinterpret_cast<SiteModel *>(siteModel_->GetManagedObject()), sitePattern, bm,
        use_tip_states);
    model_ = new_TreeLikelihoodModel("treelike", tlk, treeModel_->GetModel(),
                                     substitutionModel_->GetModel(),
                                     siteModel_->GetModel(), mbm);

    tlk->include_jacobian = include_jacobian;
    RequestGradient();
}

void TreeLikelihoodInterface::RequestGradient(
    std::vector<TreeLikelihoodGradientFlags> flags) {
    int flags_int = 0;
    for (auto flag : flags) {
        flags_int |=
            static_cast<std::underlying_type<TreeLikelihoodGradientFlags>::type>(flag);
    }
    gradientLength_ = TreeLikelihood_initialize_gradient(model_, flags_int);
    if (branchModel_ == nullptr) {
        gradientLength_ -= 2;
    }
}

void TreeLikelihoodInterface::Gradient(double *gradient) {
    double *gradient_physher = TreeLikelihood_gradient(model_);
    size_t i = 0;
    size_t j = 0;
    size_t gradientLength = gradientLength_;
    if (branchModel_ == nullptr) {
        gradientLength += 2;
        for (; i < treeModel_->GetNodeCount() - 2; i++) {
            gradient[treeModel_->nodeMap_[i]] = gradient_physher[i];
        }
        i += 2;
        j = treeModel_->GetNodeCount() - 2;
    }
    for (; i < gradientLength; i++, j++) {
        gradient[j] = gradient_physher[i];
    }
}

void TreeLikelihoodInterface::EnableSSE(bool flag) {
    SingleTreeLikelihood *tlk = reinterpret_cast<SingleTreeLikelihood *>(model_->obj);
    SingleTreeLikelihood_enable_SSE(tlk, flag);
}

CTMCScaleModelInterface::CTMCScaleModelInterface(const std::vector<double> rates,
                                                 TreeModelInterface *treeModel)
    : treeModel_(treeModel) {
    Parameters *rates_param = new_Parameters(rates.size());
    for (size_t i = 0; i < rates.size(); i++) {
        Parameters_move(rates_param, new_Parameter("ctmcscale.x", rates[i],
                                                   new_Constraint(0, INFINITY)));
    }
    ctmcScale_ = new_CTMCScale_with_parameters(
        rates_param, reinterpret_cast<Tree *>(treeModel_->GetManagedObject()));
    model_ = new_CTMCScaleModel("ctmcscale", ctmcScale_, treeModel_->GetModel());
    RequestGradient();
    parameterCount_ = 0;
}

void CTMCScaleModelInterface::RequestGradient(std::vector<GradientFlags> flags) {
    int flags_int = 0;
    for (auto flag : flags) {
        flags_int |= static_cast<std::underlying_type<GradientFlags>::type>(flag);
    }
    gradientLength_ = DistributionModel_initialize_gradient(model_, flags_int);
}

void CTMCScaleModelInterface::Gradient(double *gradient) {
    Tree_update_heights(reinterpret_cast<Tree *>(treeModel_->GetManagedObject()));
    double *gradient_physher = DistributionModel_gradient(model_);
    memcpy(gradient, gradient_physher, sizeof(double) * gradientLength_);
}

void CTMCScaleModelInterface::SetParameters(const double *parameters) {
    Parameters_set_values(ctmcScale_->x, parameters);
}

void CoalescentModelInterface::RequestGradient(std::vector<GradientFlags> flags) {
    int flags_int = 0;
    for (auto flag : flags) {
        flags_int |= static_cast<std::underlying_type<GradientFlags>::type>(flag);
    }
    gradientLength_ = Coalescent_initialize_gradient(model_, flags_int);
}

void CoalescentModelInterface::Gradient(double *gradient) {
    double *gradient_physher = Coalescent_gradient(model_);
    memcpy(gradient, gradient_physher, sizeof(double) * gradientLength_);
}

void CoalescentModelInterface::SetParameters(const double *parameters) {
    Parameters_set_values(coalescent_->p, parameters);
}
void CoalescentModelInterface::GetParameters(double *parameters) {
    Parameters_store_value(coalescent_->p, parameters);
}

ConstantCoalescentModelInterface::ConstantCoalescentModelInterface(
    double theta, TreeModelInterface *treeModel)
    : CoalescentModelInterface(treeModel) {
    Parameter *theta_param =
        new_Parameter("constant.theta", theta, new_Constraint(0, INFINITY));
    coalescent_ = new_ConstantCoalescent(
        reinterpret_cast<Tree *>(treeModel->GetManagedObject()), theta_param);
    free_Parameter(theta_param);
    model_ =
        new_CoalescentModel2("constant", coalescent_, treeModel->GetModel(), nullptr);
    RequestGradient();
    parameterCount_ = 1;
}

PiecewiseConstantCoalescentInterface::PiecewiseConstantCoalescentInterface(
    const std::vector<double> &thetas, TreeModelInterface *treeModel)
    : CoalescentModelInterface(treeModel) {
    Parameters *thetas_param = new_Parameters(thetas.size());
    for (const auto &theta : thetas) {
        Parameters_move(thetas_param, new_Parameter("slyride.theta", theta,
                                                    new_Constraint(0, INFINITY)));
    }
    coalescent_ =
        new_SkyrideCoalescent(reinterpret_cast<Tree *>(treeModel->GetManagedObject()),
                              thetas_param, COALESCENT_THETA);
    free_Parameters(thetas_param);
    model_ =
        new_CoalescentModel2("slyride", coalescent_, treeModel->GetModel(), nullptr);
    RequestGradient();
    parameterCount_ = thetas.size();
}

PiecewiseConstantCoalescentGridInterface::PiecewiseConstantCoalescentGridInterface(
    const std::vector<double> &thetas, TreeModelInterface *treeModel, double cutoff)
    : CoalescentModelInterface(treeModel) {
    Parameters *thetas_param = new_Parameters(thetas.size());
    for (const auto &theta : thetas) {
        Parameters_move(thetas_param, new_Parameter("skygrid.theta", theta,
                                                    new_Constraint(0, INFINITY)));
    }
    coalescent_ =
        new_GridCoalescent(reinterpret_cast<Tree *>(treeModel->GetManagedObject()),
                           thetas_param, thetas.size(), cutoff, COALESCENT_THETA);
    free_Parameters(thetas_param);
    model_ =
        new_CoalescentModel2("skygrid", coalescent_, treeModel->GetModel(), nullptr);
    RequestGradient();
    parameterCount_ = thetas.size();
}
