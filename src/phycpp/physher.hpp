// Copyright 2016-2022 Mathieu Fourment.
// physher is free software under the GPLv2; see LICENSE file for details.

#pragma once

#include <stddef.h>
#include <optional>
#include <string>
#include <utility>
#include <vector>

extern "C" {
#include "phyc/demographicmodels.h"
#include "phyc/gradient.h"
#include "phyc/treelikelihood.h"
}

enum class GradientFlags {
    TREE_RATIO = GRADIENT_FLAG_TREE_RATIOS,
    TREE_HEIGHT = GRADIENT_FLAG_TREE_HEIGHTS,
    COALESCENT_THETA = GRADIENT_FLAG_COALESCENT_THETA
};

enum class TreeLikelihoodGradientFlags {
    TREE_HEIGHT = TREELIKELIHOOD_FLAG_TREE,
    SITE_MODEL = TREELIKELIHOOD_FLAG_SITE_MODEL,
    SUBSTITUTION_MODEL = TREELIKELIHOOD_FLAG_SUBSTITUTION_MODEL,
    BRANCH_MODEL = TREELIKELIHOOD_FLAG_BRANCH_MODEL
};

class ModelInterface {
   public:
    virtual void SetParameters(const double *parameters) = 0;
    virtual void GetParameters(double *parameters) = 0;

    Model *GetModel() { return model_; }

    void *GetManagedObject() { return model_->obj; }

    size_t parameterCount_;

   protected:
    Model *model_;
};

class CallableModelInterface : public ModelInterface {
   public:
    double LogLikelihood() { return model_->logP(model_); }

    virtual void Gradient(double *gradient) = 0;

    size_t gradientLength_;
};

class TreeModelInterface : public ModelInterface {
   public:
    TreeModelInterface(const std::string &newick, const std::vector<std::string> &taxa,
                       std::optional<const std::vector<double>> dates);

    virtual ~TreeModelInterface() = default;

    size_t GetNodeCount() { return nodeCount_; }

    size_t GetTipCount() { return tipCount_; }

    Tree *GetTree() { return tree_; }

   protected:
    void InitializeMap(const std::vector<std::string> &taxa);

   public:
    std::vector<size_t> nodeMap_;

   protected:
    Tree *tree_;
    size_t nodeCount_;
    size_t tipCount_;
};

class UnRootedTreeModelInterface : public TreeModelInterface {
   public:
    UnRootedTreeModelInterface(const std::string &newick,
                               const std::vector<std::string> &taxa);

    virtual ~UnRootedTreeModelInterface() { model_->free(model_); }

    void SetParameters(const double *parameters) override;

    void GetParameters(double *parameters) override {}
};

class ReparameterizedTimeTreeModelInterface : public TreeModelInterface {
   public:
    ReparameterizedTimeTreeModelInterface(const std::string &newick,
                                          const std::vector<std::string> &taxa,
                                          const std::vector<double> dates);

    virtual ~ReparameterizedTimeTreeModelInterface() { model_->free(model_); }

    void SetParameters(const double *parameters) override;

    void GetParameters(double *parameters) override {}

    void GetNodeHeights(double *heights);

    void GradientTransformJVP(double *gradient, const double *height_gradient);

    void GradientTransformJVP(double *gradient, const double *height_gradient,
                              const double *heights);

    void GradientTransformJacobian(double *gradient);

    double TransformJacobian();

   protected:
    Model *transformModel_;
};

class BranchModelInterface : public ModelInterface {
   public:
    virtual ~BranchModelInterface() = default;

    void SetParameters(const double *parameters) override;

    void GetParameters(double *parameters) override;

    void SetRates(const double *rates);

   protected:
    BranchModel *branchModel_;
};

class StrictClockModelInterface : public BranchModelInterface {
   public:
    StrictClockModelInterface(double rate, TreeModelInterface *treeModel);

    virtual ~StrictClockModelInterface() { model_->free(model_); }

    void SetRate(double rate);
};

class SimpleClockModelInterface : public BranchModelInterface {
   public:
    SimpleClockModelInterface(const std::vector<double> &rates,
                              TreeModelInterface *treeModel);

    virtual ~SimpleClockModelInterface() { model_->free(model_); }
};

class SubstitutionModelInterface : public ModelInterface {
   public:
    virtual ~SubstitutionModelInterface() = default;

   protected:
    Model *Initialize(const std::string &name, Parameters *rates, Model *frequencies,
                      Model *rates_model);

    SubstitutionModel *substModel_;
};

class JC69Interface : public SubstitutionModelInterface {
   public:
    JC69Interface();

    virtual ~JC69Interface() { model_->free(model_); }

    void SetParameters(const double *parameters) override {}

    void GetParameters(double *parameters) override {}
};

class HKYInterface : public SubstitutionModelInterface {
   public:
    HKYInterface(double kappa, const std::vector<double> &frequencies);

    virtual ~HKYInterface() { model_->free(model_); }

    void SetKappa(double kappa);

    void SetFrequencies(const double *frequencies);

    void SetParameters(const double *parameters) override;

    void GetParameters(double *parameters) override {}
};

class GTRInterface : public SubstitutionModelInterface {
   public:
    GTRInterface(const std::vector<double> &rates,
                 const std::vector<double> &frequencies);

    virtual ~GTRInterface() { model_->free(model_); }

    void SetRates(const double *rates);

    void SetFrequencies(const double *frequencies);

    void SetParameters(const double *parameters) override;

    void GetParameters(double *parameters) override {}
};

class SiteModelInterface : public ModelInterface {
   public:
    virtual ~SiteModelInterface() = default;

    void SetMu(double mu);

   protected:
    SiteModel *siteModel_;
};

class ConstantSiteModelInterface : public SiteModelInterface {
   public:
    explicit ConstantSiteModelInterface(std::optional<double> mu);

    virtual ~ConstantSiteModelInterface() { model_->free(model_); }

    void SetParameters(const double *parameters) override;

    void GetParameters(double *parameters) override;
};

class WeibullSiteModelInterface : public SiteModelInterface {
   public:
    WeibullSiteModelInterface(double shape, size_t categories,
                              std::optional<double> mu);

    virtual ~WeibullSiteModelInterface() { model_->free(model_); }

    void SetShape(double shape);

    void SetParameters(const double *parameters) override;

    void GetParameters(double *parameters) override;
};

class TreeLikelihoodInterface : public CallableModelInterface {
   public:
    TreeLikelihoodInterface(
        const std::vector<std::pair<std::string, std::string>> &alignment,
        TreeModelInterface *treeModel, SubstitutionModelInterface *substitutionModel,
        SiteModelInterface *siteModel,
        std::optional<BranchModelInterface *> branchModel, bool use_ambiguities = false,
        bool use_tip_states = false, bool include_jacobian = false);

    virtual ~TreeLikelihoodInterface() { model_->free(model_); }

    void RequestGradient(std::vector<TreeLikelihoodGradientFlags> flags =
                             std::vector<TreeLikelihoodGradientFlags>());

    void Gradient(double *gradient) override;

    void SetParameters(const double *parameters) override{};

    void GetParameters(double *parameters) override{};

   private:
    TreeModelInterface *treeModel_;
    SubstitutionModelInterface *substitutionModel_;
    SiteModelInterface *siteModel_;
    BranchModelInterface *branchModel_;
    SitePattern *sitePattern_;
};

class CTMCScaleModelInterface : public CallableModelInterface {
   public:
    CTMCScaleModelInterface(const std::vector<double> rates,
                            TreeModelInterface *treeModel);

    virtual ~CTMCScaleModelInterface() = default;

    void RequestGradient(
        std::vector<GradientFlags> flags = std::vector<GradientFlags>());

    void Gradient(double *gradient) override;

    void SetParameters(const double *parameters) override;

    void GetParameters(double *parameters) override {}

   protected:
    DistributionModel *ctmcScale_;

   private:
    TreeModelInterface *treeModel_;
};

class CoalescentModelInterface : public CallableModelInterface {
   protected:
    explicit CoalescentModelInterface(TreeModelInterface *treeModel)
        : treeModel_(treeModel) {}

   public:
    virtual ~CoalescentModelInterface() = default;

    void RequestGradient(
        std::vector<GradientFlags> flags = std::vector<GradientFlags>());

    void Gradient(double *gradient) override;

    void SetParameters(const double *parameters) override;

    void GetParameters(double *parameters) override;

   protected:
    Coalescent *coalescent_;

   private:
    TreeModelInterface *treeModel_;
};

class ConstantCoalescentModelInterface : public CoalescentModelInterface {
   public:
    ConstantCoalescentModelInterface(double theta, TreeModelInterface *treeModel);

    virtual ~ConstantCoalescentModelInterface() = default;
};

class PiecewiseConstantCoalescentInterface : public CoalescentModelInterface {
   public:
    PiecewiseConstantCoalescentInterface(const std::vector<double> &thetas,
                                         TreeModelInterface *treeModel);

    virtual ~PiecewiseConstantCoalescentInterface() = default;
};

class PiecewiseConstantCoalescentGridInterface : public CoalescentModelInterface {
   public:
    PiecewiseConstantCoalescentGridInterface(const std::vector<double> &thetas,
                                             TreeModelInterface *treeModel,
                                             double cutoff);

    virtual ~PiecewiseConstantCoalescentGridInterface() = default;
};
