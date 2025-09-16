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

TimeTreeModelInterface::TimeTreeModelInterface(const std::string &newick,
                                               const std::vector<std::string> &taxa,
                                               const std::vector<double> dates) {
    char **taxa_p = new char *[taxa.size()];
    for (size_t i = 0; i < taxa.size(); i++) {
        taxa_p[i] = const_cast<char *>(taxa[i].c_str());
    }
    model_ = new_TimeTreeModel_from_newick(newick.c_str(), taxa_p, dates.data());
    delete[] taxa_p;

    tree_ = reinterpret_cast<Tree *>(model_->obj);
    nodeCount_ = Tree_node_count(tree_);
    tipCount_ = Tree_tip_count(tree_);

    InitializeMap(taxa);
    parameterCount_ = tipCount_ - 1;
}

void TimeTreeModelInterface::SetParameters(const double *parameters) {
    Node **nodes = Tree_nodes(tree_);
    for (size_t i = 0; i < nodeCount_; i++) {
        if (!Node_isleaf(nodes[i])) {
            Node_set_height_quietly(nodes[i], parameters[nodes[i]->class_id]);
        }
    }
    model_->listeners->fire(model_->listeners, model_, nullptr, -1);
}

void TimeTreeModelInterface::GetParameters(double *parameters) {
    Node **nodes = Tree_nodes(tree_);
    for (size_t i = 0; i < nodeCount_; i++) {
        if (!Node_isleaf(nodes[i])) {
            parameters[nodes[i]->class_id] = Node_height(nodes[i]);
        }
    }
}

void TimeTreeModelInterface::GetNodeHeights(double *parameters) {
    this->GetParameters(parameters);
}

ReparameterizedTimeTreeModelInterface::ReparameterizedTimeTreeModelInterface(
    const std::string &newick, const std::vector<std::string> &taxa,
    const std::vector<double> dates, TreeTransformFlags transform)
    : TimeTreeModelInterface(newick, taxa, dates) {
    int int_transform =
        static_cast<std::underlying_type<TreeLikelihoodGradientFlags>::type>(transform);
    TreeModel_set_transform(model_, int_transform);
    transformModel_ = reinterpret_cast<Model *>(model_->data);
}

void ReparameterizedTimeTreeModelInterface::SetParameters(const double *parameters) {
    TreeTransform *tt = reinterpret_cast<TreeTransform *>(transformModel_->obj);
    Parameters_set_values_quietly(tt->parameters, parameters);
    transformModel_->listeners->fire(transformModel_->listeners, transformModel_, nullptr, -1);
}

void ReparameterizedTimeTreeModelInterface::GetParameters(double *parameters) {
    TreeTransform *tt = reinterpret_cast<TreeTransform *>(transformModel_->obj);
    Parameters_set_values(tt->parameters, parameters);
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
    DiscreteParameter *map = new_DiscreteParameter(nodeCount);
    Parameter* rates_parameter = new_Parameter2("clock_rates", rates.data(), rates.size(), new_Constraint(0, INFINITY));
    Tree *tree = treeModel->GetTree();
    for (size_t i = 0; i < nodeCount; i++) {
        Node *node = Tree_node(tree, i);
        if (!Node_isroot(node)) {
            size_t index = Node_isleaf(node)
                               ? node->class_id
                               : node->class_id + treeModel->GetTipCount();
            map->values[Node_id(node)] = index;
        }
    }
    parameterCount_ = nodeCount - 1;
    branchModel_ = new_DiscreteClock_with_parameters(treeModel->GetTree(), rates_parameter, map);
    model_ =
        new_BranchModel2("simple_clock", branchModel_, treeModel->GetModel(), nullptr);
}

JC69Interface::JC69Interface() {
    double values[4] = {0.25, 0.25, 0.25, 0.25};
    Parameter* frequencies = new_Parameter2("jc69_frequencies", values, 4, new_Constraint(0, 1));
    dataType_ = new NucleotideDataTypeInterface();
    substModel_ = new_JC69(frequencies);
    free_DataType(substModel_->datatype);
    substModel_->datatype = dataType_->dataType_;
    substModel_->datatype->ref_count++;
    model_ = new_SubstitutionModel2("jc69", substModel_);
    parameterCount_ = 0;
}

HKYInterface::HKYInterface(double kappa, const std::vector<double> &frequencies) {
    Parameter *kappa_parameter =
        new_Parameter("hky_kappa", kappa, new_Constraint(0, INFINITY));
    Parameter* frequencies_parameter = new_Parameter2("hky_frequencies", frequencies.data(), 4, new_Constraint(0, 1));
    dataType_ = new NucleotideDataTypeInterface();
    substModel_ = new_HKY_with_parameters(frequencies_parameter, kappa_parameter);
    free_DataType(substModel_->datatype);
    substModel_->datatype = dataType_->dataType_;
    substModel_->datatype->ref_count++;
    model_ = new_SubstitutionModel2("hky", substModel_);
    free_Parameter(kappa_parameter);
    parameterCount_ = 5;
}

void HKYInterface::SetKappa(double kappa) {
    Parameters_set_value(substModel_->rates, 0, kappa);
}

void HKYInterface::SetFrequencies(const double *frequencies) {
    Parameter_set_values(substModel_->simplex, frequencies);
}

void HKYInterface::SetParameters(const double *parameters) {
    SetKappa(parameters[0]);
    Parameter_set_values(substModel_->simplex, parameters + 1);
}

GTRInterface::GTRInterface(const std::vector<double> &rates,
                           const std::vector<double> &frequencies) {
    Parameters *rates_parameters = nullptr;
    Parameter *rates_simplex = nullptr;
    // simplex
    if (rates.size() == 6) {
        rates_simplex = new_Parameter2("gtr_rate", rates.data(), 6, new_Constraint(0, 1));
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
    Parameter* frequencies_parameter = new_Parameter2("gtr_frequencies", frequencies.data(), 4, new_Constraint(0, 1));

    dataType_ = new NucleotideDataTypeInterface();
    const char *assigments[5] = {"AC", "AG", "AT", "CG", "CT"};

    substModel_ =
        SubstitutionModel_factory("gtr", dataType_->dataType_,
                                  frequencies_parameter,
                                  rates_simplex, rates_parameters, assigments);
    free_DataType(substModel_->datatype);
    substModel_->datatype = dataType_->dataType_;
    substModel_->datatype->ref_count++;
    model_ = new_SubstitutionModel2("gtr", substModel_);
    free_Parameters(rates_parameters);
}

void GTRInterface::SetRates(const double *rates) {
    if (substModel_->rates == nullptr) {
        Parameter_set_values(substModel_->rates_simplex, rates);
    } else {
        Parameters_set_values(substModel_->rates, rates);
    }
}
void GTRInterface::SetFrequencies(const double *frequencies) {
    Parameter_set_values(substModel_->simplex, frequencies);
}
void GTRInterface::SetParameters(const double *parameters) {
    SetRates(parameters);
    size_t offset = (substModel_->rates == nullptr ? 6 : 5);
    SetFrequencies(parameters + offset);
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
    Parameter* frequencies_parameter = new_Parameter2("general_frequencies", frequencies.data(), frequencies.size(), new_Constraint(0, 1));

    DiscreteParameter *dp =
        new_DiscreteParameter_with_values(structure.data(), structure.size());
    Model *mdp = new_DiscreteParameterModel("structure", dp);

    substModel_ = new_GeneralModel_with_parameters(
        dataType->dataType_, reinterpret_cast<DiscreteParameter *>(mdp->obj),
        rate_parameters, frequencies_parameter, -1, normalize);
    model_ = new_SubstitutionModel3("substmodel", substModel_, mdp);

    free_Parameters(rate_parameters);
    dataType_ = dataType;
}

void GeneralSubstitutionModelInterface::SetRates(const double *rates) {
        Parameters_set_values(substModel_->rates, rates);
}

void GeneralSubstitutionModelInterface::SetFrequencies(const double *frequencies) {
    Parameter_set_values(substModel_->simplex, frequencies);
}

void GeneralSubstitutionModelInterface::SetParameters(const double *parameters) {
    Parameters_set_values(substModel_->rates, parameters);
    Parameter_set_values(substModel_->simplex, parameters + Parameters_count(substModel_->rates));
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
    model_ = new_SiteModel2("sitemodel", siteModel_);
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
void ConstantSiteModelInterface::GetRates(double *rates) {
    rates[0] = (siteModel_->mu != nullptr ? Parameter_value(siteModel_->mu) : 1.0);
}

void ConstantSiteModelInterface::GetProportions(double *proportions) {
    proportions[0] = 1.0;
}

InvariantSiteModelInterface::InvariantSiteModelInterface(double proportionInvariant,
                                                         std::optional<double> mu) {
    Parameter* prop_parameter = new_Parameter("general_frequencies", proportionInvariant, new_Constraint(0, 1));

    siteModel_ =
        new_SiteModel_with_parameters(nullptr, prop_parameter, 1, DISTRIBUTION_DISCRETE,
                                      true, QUADRATURE_QUANTILE_MEDIAN);
    if (mu.has_value()) {
        Parameter *mu_param = new_Parameter("mu", *mu, new_Constraint(0, INFINITY));
        SiteModel_set_mu(siteModel_, mu_param);
        free_Parameter(mu_param);
    }

    model_ = new_SiteModel2("sitemodel", siteModel_);
    parameterCount_ = 1 + siteModel_->mu != nullptr;
}

void InvariantSiteModelInterface::SetProportionInvariant(double value) {
    Parameter_set_value(siteModel_->proportions, value);
}

void InvariantSiteModelInterface::SetParameters(const double *parameters) {
    size_t shift = 0;
    if (siteModel_->proportions != nullptr) {
        Parameter_set_value(siteModel_->proportions, *parameters);
        shift++;
    }

    if (siteModel_->mu != nullptr) {
        Parameter_set_value(siteModel_->mu, parameters[shift]);
    }
}

void InvariantSiteModelInterface::GetParameters(double *parameters) {
    size_t shift = 0;
    if (siteModel_->proportions != nullptr) {
        parameters[0] = Parameter_value(siteModel_->proportions);
        shift++;
    }
    if (siteModel_->mu != nullptr) {
        parameters[shift] = Parameter_value(siteModel_->mu);
    }
}

void InvariantSiteModelInterface::GetRates(double *rates) {
    siteModel_->update(siteModel_);
    memcpy(rates, siteModel_->cat_rates, sizeof(double) * siteModel_->cat_count);
}

void InvariantSiteModelInterface::GetProportions(double *proportions) {
    siteModel_->update(siteModel_);
    memcpy(proportions, siteModel_->cat_proportions,
           sizeof(double) * siteModel_->cat_count);
}

DiscretizedSiteModelInterface::DiscretizedSiteModelInterface(
    distribution_t distribution, double shape, size_t categories,
    std::optional<double> proportionInvariant, std::optional<double> mu) {
    Parameter *shape_parameter =
        new_Parameter("sitemodel.parameter", shape, new_Constraint(0, INFINITY));
    Parameters *params = new_Parameters(1);
    Parameters_move(params, shape_parameter);

    Parameter *propInvSimplex = nullptr;
    if (proportionInvariant.has_value()) {
        double values[2] = {proportionInvariant.value(), 1.0 - proportionInvariant.value()};
        propInvSimplex = new_Parameter2("proportionInvariant", values, 2, new_Constraint(0, 1));
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

    model_ = new_SiteModel2("sitemodel", siteModel_);
    free_Parameters(params);
    parameterCount_ = Parameters_count(siteModel_->rates) + siteModel_->mu != nullptr;
    if (propInvSimplex != nullptr) {
        parameterCount_++;
    }
    categoryCount_ = siteModel_->cat_count;
}

void DiscretizedSiteModelInterface::SetParameter(double parameter) {
    Parameters_set_value(siteModel_->rates, 0, parameter);
}

void DiscretizedSiteModelInterface::SetProportionInvariant(double value) {
    double values[2] = {value, 1.0 - value};
    Parameter_set_values(siteModel_->proportions, values);
}

void DiscretizedSiteModelInterface::SetParameters(const double *parameters) {
    size_t shift = 0;
    if (Parameters_count(siteModel_->rates) > 0) {
        Parameters_set_value(siteModel_->rates, 0, parameters[0]);
        shift++;
    }
    if (siteModel_->proportions != nullptr) {
        double values[2] = {parameters[shift], 1.0 - parameters[shift]};
        Parameter_set_values(siteModel_->proportions, values);
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
        parameters[shift] = Parameter_value(siteModel_->proportions);
        shift++;
    }
    if (siteModel_->mu != nullptr) {
        parameters[shift] = Parameter_value(siteModel_->mu);
        shift++;
    }
}

void DiscretizedSiteModelInterface::GetRates(double *rates) {
    siteModel_->update(siteModel_);
    memcpy(rates, siteModel_->cat_rates, sizeof(double) * siteModel_->cat_count);
}

void DiscretizedSiteModelInterface::GetProportions(double *proportions) {
    siteModel_->update(siteModel_);
    memcpy(proportions, siteModel_->cat_proportions,
           sizeof(double) * siteModel_->cat_count);
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
}

void TreeLikelihoodInterface::Gradient(double *gradient, std::vector<GradientFlags> flags) {
    int treelikelihoodFlags = 0;
    for(auto flag : flags) {
        switch(flag){
            case GradientFlags::BRANCH_MODEL_RATE:
                treelikelihoodFlags |= TREELIKELIHOOD_FLAG_BRANCH_MODEL;
                break;
            case GradientFlags::SITE_MODEL_MU:
            case GradientFlags::SITE_MODEL_PARAMETER:
            case GradientFlags::SITE_MODEL_PINV:
                treelikelihoodFlags |= TREELIKELIHOOD_FLAG_SITE_MODEL;
                break;
            case GradientFlags::SUBSTITUTION_MODEL_FREQUENCIES:
                treelikelihoodFlags |= TREELIKELIHOOD_FLAG_SUBSTITUTION_MODEL_FREQUENCIES;
                break;
            case GradientFlags::SUBSTITUTION_MODEL_RATES:
                treelikelihoodFlags |= TREELIKELIHOOD_FLAG_SUBSTITUTION_MODEL_RATES;
                break;
            case GradientFlags::TREE_MODEL_BRANCH_LENGTHS:
            case GradientFlags::TREE_MODEL_HEIGHT:
            case GradientFlags::TREE_MODEL_RATIO:
                treelikelihoodFlags |= TREELIKELIHOOD_FLAG_TREE_MODEL;
                break;
            default:
                std::cerr << "TreeLikelihoodInterface: unknown gradient flag " << static_cast<std::underlying_type<GradientFlags>::type>(flag) << std::endl;
        }
    }
    TreeLikelihood_gradient(model_, treelikelihoodFlags, gradient);
}

void TreeLikelihoodInterface::EnableSSE(bool flag) {
    SingleTreeLikelihood *tlk = reinterpret_cast<SingleTreeLikelihood *>(model_->obj);
    SingleTreeLikelihood_enable_SSE(tlk, flag);
}

int CallableModelInterface::ConvertFlags(std::vector<GradientFlags> flags){
    int flagInts = 0;
    for (auto flag : flags) {
        flagInts |= static_cast<std::underlying_type<GradientFlags>::type>(flag);
    }
    return flagInts;
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
    parameterCount_ = 0;
}

void CTMCScaleModelInterface::Gradient(double *gradient, std::vector<GradientFlags> flags) {
    int flagInts = this->ConvertFlags(flags);
    CTMCModel_gradient(model_, flagInts, gradient);
}

void CTMCScaleModelInterface::SetParameters(const double *parameters) {
    Parameters_set_values(ctmcScale_->x, parameters);
}

void CTMCScaleModelInterface::GetParameters(double *parameters) {
    Parameters_store_value(ctmcScale_->x, parameters);
}

void CoalescentModelInterface::Gradient(double *gradient, std::vector<GradientFlags> flags) {
   int flagInts = this->ConvertFlags(flags);
    Coalescent_gradient(model_, flagInts, gradient);
}

void CoalescentModelInterface::SetParameters(const double *parameters) {
    Parameters_set_values(coalescent_->p, parameters);
}
void CoalescentModelInterface::GetParameters(double *parameters) {
    Parameters_store_value(coalescent_->p, parameters);
}

ConstantCoalescentModelInterface::ConstantCoalescentModelInterface(
    double theta, TimeTreeModelInterface *treeModel)
    : CoalescentModelInterface(treeModel) {
    Parameter *theta_param =
        new_Parameter("constant.theta", theta, new_Constraint(0, INFINITY));
    coalescent_ = new_ConstantCoalescent(
        reinterpret_cast<Tree *>(treeModel->GetManagedObject()), theta_param);
    free_Parameter(theta_param);
    model_ =
        new_CoalescentModel2("constant", coalescent_, treeModel->GetModel(), nullptr);
    parameterCount_ = 1;
}

PiecewiseConstantCoalescentInterface::PiecewiseConstantCoalescentInterface(
    const std::vector<double> &thetas, TimeTreeModelInterface *treeModel)
    : CoalescentModelInterface(treeModel) {
    Parameter *thetas_param = new_Parameter2("skyride.theta", thetas.data(), thetas.size(),
                                                    new_Constraint(0, INFINITY));
    coalescent_ =
        new_SkyrideCoalescent(reinterpret_cast<Tree *>(treeModel->GetManagedObject()),
                              thetas_param);
    // free_Parameters(thetas_param);
    model_ =
        new_CoalescentModel2("slyride", coalescent_, treeModel->GetModel(), nullptr);
    parameterCount_ = thetas.size();
}

PiecewiseConstantCoalescentGridInterface::PiecewiseConstantCoalescentGridInterface(
    const std::vector<double> &thetas, TimeTreeModelInterface *treeModel, double cutoff)
    : CoalescentModelInterface(treeModel) {
        Parameter *thetas_param = new_Parameter2("skygrid.theta", thetas.data(), thetas.size(),
                                                    new_Constraint(0, INFINITY));
    coalescent_ =
        new_GridCoalescent(reinterpret_cast<Tree *>(treeModel->GetManagedObject()),
                           thetas_param, thetas.size(), cutoff);
    // free_Parameters(thetas_param);
    model_ =
        new_CoalescentModel2("skygrid", coalescent_, treeModel->GetModel(), nullptr);
    parameterCount_ = thetas.size();
}

PiecewiseLinearCoalescentGridInterface::PiecewiseLinearCoalescentGridInterface(
    const std::vector<double> &thetas, TimeTreeModelInterface *treeModel, double cutoff)
    : CoalescentModelInterface(treeModel) {
        Parameter *thetas_param = new_Parameter2("piecewise-linear.theta", thetas.data(), thetas.size(),
                                                    new_Constraint(0, INFINITY));
    coalescent_ = new_PiecewiseLinearGridCoalescent(
        reinterpret_cast<Tree *>(treeModel->GetManagedObject()), thetas_param,
        thetas.size(), cutoff);
    // free_Parameters(thetas_param);
    model_ = new_CoalescentModel2("piecewise-linear", coalescent_,
                                  treeModel->GetModel(), nullptr);
    parameterCount_ = thetas.size();
}
